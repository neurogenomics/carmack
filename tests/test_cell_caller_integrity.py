"""Small independently specified fixtures for exported matrix identity and CIGAR overlap."""

import pysam
import pytest
from scipy.io import mmread

from carmack.cell_caller.cell_caller import CellCaller
from carmack.cell_caller.peak_barcode_matrix import PeakBarcodeMatrix


def make_caller(tmp_path, bed_text, reads):
    """Write sorted alignments; each tuple is (start, CIGAR, query length)."""
    bed = tmp_path / "peaks.bed"
    bed.write_text(bed_text)
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i, (start, cigar, length) in enumerate(reads):
            read = pysam.AlignedSegment(out.header)
            read.query_name = f"read{i}"
            read.query_sequence = "A" * length
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = cigar
            read.set_tag("CB", "cell")
            out.write(read)
    pysam.index(str(bam))
    return CellCaller(str(bed), str(bam))


def test_export_peak_coordinates_follow_matrix_row_order(tmp_path):
    caller = make_caller(
        tmp_path,
        "chr1\t100\t120\tz_peak\t0\t.\nchr1\t300\t320\ta_peak\t0\t.\n",
        [(100, "20M", 20), (300, "20M", 20), (300, "20M", 20)],
    )
    caller.compute_matrix()
    caller.export(str(tmp_path), "out", force_n=0)
    counts = mmread(str(tmp_path / "out_matrix.mtx")).toarray()[:, 0]
    coordinates = [
        tuple(line.split("\t")) for line in (tmp_path / "out_peaks.bed").read_text().splitlines()
    ]
    assert dict(zip(coordinates, counts)) == {("chr1", "100", "120"): 1, ("chr1", "300", "320"): 2}


def test_repeated_interval_is_exported_once(tmp_path):
    caller = make_caller(
        tmp_path,
        "chr1\t100\t120\tpeak\t0\t.\nchr1\t100\t120\tpeak\t0\t.\n",
        [(100, "20M", 20)],
    )
    caller.compute_matrix()
    caller.export(str(tmp_path), "out", force_n=0)
    coordinates = (tmp_path / "out_peaks.bed").read_text().splitlines()
    matrix = mmread(str(tmp_path / "out_matrix.mtx"))
    assert len(coordinates) == matrix.shape[0] == 1
    assert matrix.toarray().tolist() == [[1]]


def test_distinct_intervals_cannot_share_peak_name(tmp_path):
    caller = make_caller(
        tmp_path,
        "chr1\t100\t120\tpeak\t0\t.\nchr1\t300\t320\tpeak\t0\t.\n",
        [(100, "20M", 20), (300, "20M", 20)],
    )
    with pytest.raises(ValueError, match="[Dd]uplicate peak"):
        caller.compute_matrix()


@pytest.mark.parametrize("cigar", ["10M100N10M", "10M100D10M"])
def test_unaligned_reference_gap_does_not_overlap_peak(tmp_path, cigar):
    caller = make_caller(tmp_path, "chr1\t130\t140\tgap\t0\t.\n", [(100, cigar, 20)])
    caller.compute_matrix()
    assert caller.matrix.sum_barcodes() == [0]


def test_peak_matrix_accepts_one_shot_iterables():
    matrix = PeakBarcodeMatrix((x for x in ["b", "a"]), (x for x in ["cell"]))
    matrix.increment_index("b", "cell")
    assert matrix.get_value("b", "cell") == 1
    assert matrix.matrix.shape == (2, 1)


def test_flat_rank_curve_requires_explicit_cell_count(tmp_path):
    caller = make_caller(tmp_path, "chr1\t100\t120\tpeak\t0\t.\n", [(100, "20M", 20)])
    caller.compute_matrix()
    with pytest.raises(ValueError, match="force_n"):
        caller.export(str(tmp_path), "out")
    assert not (tmp_path / "out_matrix.mtx").exists()


def test_empty_matrix_can_be_exported_with_force_all(tmp_path):
    caller = make_caller(tmp_path, "chr1\t100\t120\tpeak\t0\t.\n", [])
    caller.compute_matrix()
    caller.export(str(tmp_path), "out", force_n=0)
    assert mmread(str(tmp_path / "out_matrix.mtx")).shape == (0, 0)
    assert (tmp_path / "out_peaks.bed").read_text() == ""
    assert (tmp_path / "out_barcodes.tsv").read_text() == ""
