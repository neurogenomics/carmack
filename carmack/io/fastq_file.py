from functools import cached_property
from itertools import islice

from .gzip_file import GzipFile
from .subprocess_stream import SubprocessStream

GZIP_SUFFIX = ".gz"
LZ4_SUFFIX = ".lz4"
LINES_PER_READ = 4
LINES_PER_PAIRED_READ = 8


class FastqFile(GzipFile):
    """
    Class that can read/write fastq files in raw or gz format
    """

    def __init__(self, filename: str, paired_end: bool = False):
        """
        Initialise the FastqFile object
        """
        self.filename = filename
        self.paired_end = paired_end
        super().__init__(filename)

    @cached_property
    def reads_count(self) -> int:
        """
        Count the total number of reads in the FASTQ file.

        Returns:
            int: Total number of reads in the file.
        """
        line_count = 0

        if self.compressor is not None:
            stream = SubprocessStream([self.compressor, "-c", "-d", self.filename], mode="r")
        else:
            stream = open(self.filename, "r")

        with stream as fastq_file:
            for _ in fastq_file:
                line_count += 1

        lines_per_record = LINES_PER_PAIRED_READ if self.paired_end else LINES_PER_READ
        if line_count % lines_per_record:
            raise ValueError(
                f"Invalid FASTQ {self.filename!r}: incomplete record or interleaved pair "
                f"({line_count} lines; expected a multiple of {lines_per_record})."
            )
        return line_count // lines_per_record

    def open_read_iterator(self, as_string: bool = False):
        """Yield complete four-line FASTQ records, rejecting malformed or truncated input.

        ``paired_end`` means interleaved records in one file. With ``as_string=True``
        every field is text for both plain and compressed inputs. The legacy default
        retains the underlying stream's type (text for plain, bytes for compressed).
        """
        if self.compressor is not None:
            stream = SubprocessStream([self.compressor, "-c", "-d", self.filename], mode="r")
        else:
            stream = open(self.filename, "r")

        pending_mate = None
        partial_error = None
        record_number = 0
        with stream as fastq_file:
            for header in fastq_file:
                record_number += 1
                record = [header, *islice(fastq_file, LINES_PER_READ - 1)]
                context = f"Invalid FASTQ {self.filename!r} at record {record_number}"
                if len(record) != LINES_PER_READ:
                    partial_error = ValueError(f"{context}: incomplete four-line record.")
                    break
                newline = b"\r\n" if isinstance(header, bytes) else "\r\n"
                header, seq, separator, qual = (line.rstrip(newline) for line in record)
                at = b"@" if isinstance(header, bytes) else "@"
                plus = b"+" if isinstance(header, bytes) else "+"
                if not header.startswith(at) or len(header) == 1:
                    raise ValueError(f"{context}: expected a non-empty @ read header.")
                if not separator.startswith(plus):
                    raise ValueError(f"{context}: expected a + separator.")
                if len(seq) != len(qual):
                    raise ValueError(
                        f"{context}: sequence and quality lengths differ "
                        f"({len(seq)} != {len(qual)})."
                    )
                fields = (header[1:], seq, qual)
                if as_string and isinstance(header, bytes):
                    fields = tuple(field.decode("UTF-8") for field in fields)
                if not self.paired_end:
                    yield fields
                elif pending_mate is None:
                    pending_mate = fields
                else:
                    yield (*pending_mate, *fields)
                    pending_mate = None

        # Check the decompressor's exit status first: a broken gzip stream should retain
        # its more informative compressor failure, including when it cut a record short.
        if partial_error is not None:
            raise partial_error
        if pending_mate is not None:
            raise ValueError(
                f"Invalid FASTQ {self.filename!r} at record {record_number}: "
                "interleaved input ends without the second mate."
            )

    @staticmethod
    def write_read(file_stream, name: str, seq: str, qual: str) -> None:
        """Writes a single read to a fastq file"""
        file_stream.write(("@" + name + "\n").encode("UTF-8"))
        file_stream.write((seq + "\n").encode("UTF-8"))
        file_stream.write(("+\n").encode("UTF-8"))
        file_stream.write((qual + "\n").encode("UTF-8"))
