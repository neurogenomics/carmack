"""FASTQ framing errors must not be accepted as a successful shorter run."""

import gzip

import pytest

from carmack.io.fastq_file import FastqFile


def write_fastq(tmp_path, text, compressed):
    path = tmp_path / ("reads.fastq.gz" if compressed else "reads.fastq")
    path.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
    return str(path)


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize(
    "text",
    [
        "@ok\nAC\n+\nII\n@truncated\nAC\n+\n",
        "not-a-header\nAC\n+\nII\n",
        "@read\nAC\nnot-plus\nII\n",
        "@read\nAC\n+\nI\n",
    ],
)
def test_invalid_fastq_fails_with_filename_and_record(tmp_path, text, compressed):
    path = write_fastq(tmp_path, text, compressed)
    with pytest.raises(ValueError, match=r"FASTQ.*record") as error:
        list(FastqFile(path).open_read_iterator(as_string=True))
    assert path in str(error.value)


@pytest.mark.parametrize("compressed", [False, True])
def test_read_count_rejects_partial_records(tmp_path, compressed):
    path = write_fastq(tmp_path, "@read\nAC\n+\n", compressed)
    with pytest.raises(ValueError, match="FASTQ"):
        _ = FastqFile(path).reads_count


@pytest.mark.parametrize("compressed", [False, True])
def test_interleaved_pair_decodes_plain_and_compressed(tmp_path, compressed):
    path = write_fastq(tmp_path, "@read/1\nAC\n+\nII\n@read/2\nGT\n+\nJJ\n", compressed)
    assert list(FastqFile(path, paired_end=True).open_read_iterator(as_string=True)) == [
        ("read/1", "AC", "II", "read/2", "GT", "JJ")
    ]


@pytest.mark.parametrize("compressed", [False, True])
def test_interleaved_input_cannot_end_with_unpaired_record(tmp_path, compressed):
    path = write_fastq(tmp_path, "@read/1\nAC\n+\nII\n", compressed)
    with pytest.raises(ValueError, match="FASTQ"):
        list(FastqFile(path, paired_end=True).open_read_iterator(as_string=True))
