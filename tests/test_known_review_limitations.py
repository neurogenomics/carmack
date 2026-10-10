"""Executable open review findings; xfail means unresolved, never biological validation."""

from unittest.mock import MagicMock

import pysam
import pytest

from carmack.io import log_subprocess
from carmack.linear_dedup.linear_dedup import LinearDedup
from carmack.tag_dedup.tag_dedup import TagDedup


def write_pairs(path, records):
    """Records specify name, mate number and barcode in coordinate-tie order."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, mate, barcode, score in records:
            read = pysam.AlignedSegment(out.header)
            read.query_name = name
            read.query_sequence = "A" * 50
            read.flag = 99 if mate == 1 else 147
            read.reference_id = read.next_reference_id = 0
            read.reference_start = 100 if mate == 1 else 250
            read.next_reference_start = 250 if mate == 1 else 100
            read.template_length = 200 if mate == 1 else -200
            read.cigarstring = "50M"
            read.set_tag("AS", score)
            if barcode is not None:
                read.set_tag("CB", barcode)
            out.write(read)
    pysam.index(str(path))


@pytest.mark.xfail(strict=True, reason="Open review: legacy dedup chooses mates independently")
def test_legacy_dedup_preserves_a_complete_winning_pair(tmp_path):
    bam = tmp_path / "reads.bam"
    write_pairs(bam, [("a", 1, None, 1), ("b", 1, None, 1), ("b", 2, None, 1), ("a", 2, None, 1)])
    csv = tmp_path / "barcodes.csv"
    csv.write_text("a,ACGT\nb,ACGT\n")
    TagDedup(str(bam), str(bam) + ".bai", str(csv)).tag_dedup_reads(True, str(tmp_path), "out")
    with pysam.AlignmentFile(str(tmp_path / "out.dedup.tagged.bam"), "rb") as out:
        records = [(read.query_name, read.is_read1) for read in out]
    assert records == [("a", True), ("a", False)]


@pytest.mark.xfail(strict=True, reason="Open review: linear-dedup requires globally unique QNAMEs")
def test_linear_dedup_does_not_resurrect_a_loser_with_a_shared_name(tmp_path):
    bam = tmp_path / "reads.bam"
    write_pairs(
        bam,
        [
            ("shared", 1, "A", 20),
            ("shared", 1, "B", 10),
            ("winner", 1, "B", 30),
            ("shared", 2, "A", 20),
            ("shared", 2, "B", 10),
            ("winner", 2, "B", 30),
        ],
    )
    stats = LinearDedup(str(bam), str(bam) + ".bai").linear_dedup_reads(str(tmp_path), "out")
    with pysam.AlignmentFile(str(tmp_path / "out.linear_dedup.bam"), "rb") as out:
        records = [(read.query_name, read.get_tag("CB")) for read in out]
    assert len(records) == 2 * stats.pairs_kept == 4
    assert ("shared", "B") not in records


@pytest.mark.xfail(strict=True, reason="Open review: Linux parent-death callback runs in parent")
def test_linux_subprocess_setup_defers_prctl_until_child(monkeypatch):
    libc = MagicMock()
    libc.prctl.return_value = 0
    monkeypatch.setattr(log_subprocess, "LIBC", libc)
    monkeypatch.setattr(log_subprocess.sys, "platform", "linux")
    wrapper = log_subprocess.LogSubprocess()
    libc.prctl.assert_not_called()
    assert callable(wrapper.pdeathsig)
