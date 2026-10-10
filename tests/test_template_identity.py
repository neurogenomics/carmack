"""Template-level regression fixtures, independent of biological validation."""

from collections import Counter

import pysam
import pytest

from carmack.bam_templates import PrimaryPairValidator, template_identity
from carmack.linear_dedup.linear_dedup import LinearDedup
from carmack.tag_dedup.tag_dedup import TagDedup


def pair(name, cell="A", group=None, start=100, end=250, score=20, chrom=0):
    result = []
    for mate in (1, 2):
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "A" * 50
        read.flag = 99 if mate == 1 else 147
        read.reference_id = read.next_reference_id = chrom
        read.reference_start = start if mate == 1 else end
        read.next_reference_start = end if mate == 1 else start
        read.template_length = (end + 50 - start) * (1 if mate == 1 else -1)
        read.cigarstring = "50M"
        read.set_tag("CB", cell)
        read.set_tag("AS", score)
        if group is not None:
            read.set_tag("RG", group)
        result.append(read)
    return result


def write_input(tmp_path, records, kind):
    path = tmp_path / "input.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}],
        "RG": [{"ID": "lane1"}, {"ID": "lane2"}],
    }
    barcodes = {}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for read in sorted(
            records,
            key=lambda r: (r.reference_id if r.reference_id >= 0 else 2**31, r.reference_start),
        ):
            if kind == "legacy":
                barcodes[read.query_name] = read.get_tag("CB")
                read.set_tag("CB", None)
            out.write(read)
    pysam.index(str(path))
    csv = tmp_path / "barcodes.csv"
    csv.write_text("".join(f"{name},{barcode}\n" for name, barcode in barcodes.items()))
    return path, csv


def run(tmp_path, records, kind="linear"):
    path, csv = write_input(tmp_path, records, kind)
    if kind == "linear":
        stats = LinearDedup(str(path), str(path) + ".bai").linear_dedup_reads(str(tmp_path), "out")
        output = tmp_path / "out.linear_dedup.bam"
    else:
        TagDedup(str(path), str(path) + ".bai", str(csv)).tag_dedup_reads(
            True, str(tmp_path), "out"
        )
        stats = None
        output = tmp_path / "out.dedup.tagged.bam"
    with pysam.AlignmentFile(str(output), "rb") as bam:
        return list(bam.fetch(until_eof=True)), stats


@pytest.mark.parametrize("namespace", ["cell", "readgroup", "missing_readgroup"])
def test_linear_name_reuse_cannot_resurrect_loser(tmp_path, namespace):
    if namespace == "cell":
        kept = pair("shared", cell="A")
        lost = pair("shared", cell="B", score=10)
        winner = pair("winner", cell="B", score=30)
    else:
        kept = pair("shared", group="lane1" if namespace == "readgroup" else None)
        lost = pair("shared", group="lane2", start=500, end=650, score=10)
        winner = pair("winner", group="lane2", start=500, end=650, score=30)
    output, stats = run(tmp_path, kept + lost + winner)
    assert stats.pairs_kept() == 2
    assert len(output) == 4
    shared = [r for r in output if r.query_name == "shared"]
    assert len(shared) == 2
    assert {r.reference_start for r in shared} == {100, 250}
    assert {r.get_tag("CB") for r in shared} == {"A"}
    if namespace != "cell":
        assert {r.get_tag("RG") if r.has_tag("RG") else None for r in shared} == {
            "lane1" if namespace == "readgroup" else None
        }


@pytest.mark.parametrize("kind", ["linear", "legacy"])
@pytest.mark.parametrize(
    "defect",
    [
        "repeated_mate",
        "missing_mate",
        "mate_position",
        "mate_strand",
        "both_mate_flags",
        "inconsistent_group",
        "different_contig_collision",
    ],
)
def test_ambiguous_or_incomplete_primary_pairs_fail_before_publication(tmp_path, kind, defect):
    records = pair("shared")
    if defect == "repeated_mate":
        records += pair("shared", score=10)
    elif defect == "missing_mate":
        records.pop()
    elif defect == "mate_position":
        records[1].next_reference_start += 1
    elif defect == "mate_strand":
        records[1].mate_is_reverse = not records[1].mate_is_reverse
    elif defect == "both_mate_flags":
        records[1].is_read1 = True
    elif defect == "inconsistent_group":
        records[1].set_tag("RG", "lane2")
    else:
        records += pair("shared", chrom=1)
    with pytest.raises(ValueError, match="template|Template|primary|Primary"):
        run(tmp_path, records, kind)
    assert not list(tmp_path.glob("out*.bam"))
    assert not list(tmp_path.glob("out*stats*"))


def test_linear_requires_barcode_on_both_mates(tmp_path):
    records = pair("missing_cell")
    records[1].set_tag("CB", None)
    with pytest.raises(ValueError, match="barcode tag"):
        run(tmp_path, records)
    assert not (tmp_path / "out.linear_dedup.bam").exists()


@pytest.mark.parametrize("kind", ["linear", "legacy"])
def test_readgroup_scopes_identity_but_does_not_split_duplicate_groups(tmp_path, kind):
    first = pair("same_name", group="lane1", score=30)
    second = pair("same_name", group="lane2", score=10)
    # Reverse R2 tie order is the original pair-breaking trigger.
    output, _ = run(tmp_path, [first[0], second[0], second[1], first[1]], kind)
    assert len(output) == 2
    assert {r.get_tag("RG") for r in output} == {"lane1"}
    assert {r.is_read1 for r in output} == {True, False}


@pytest.mark.parametrize("kind", ["linear", "legacy"])
def test_non_primary_records_cannot_select_a_winner(tmp_path, kind):
    first, second = pair("winner", score=30), pair("loser", score=10)
    secondary = pair("loser", start=500, end=650, score=100)[0]
    secondary.is_secondary = True
    supplementary = pair("winner", start=700, end=850, score=100)[0]
    supplementary.is_supplementary = True
    output, _ = run(tmp_path, first + second + [secondary, supplementary], kind)
    assert {r.query_name for r in output} == {"winner"}
    if kind == "linear":
        assert len(output) == 2
        assert not any(r.is_secondary or r.is_supplementary for r in output)
    else:
        assert len(output) == 3
        assert any(r.is_supplementary for r in output)


@pytest.mark.parametrize("kind", ["linear", "legacy"])
def test_unmapped_and_unpaired_records_are_not_collapsed_or_resurrected(tmp_path, kind):
    records = pair("mapped")
    for name in ("unmapped1", "unmapped2"):
        unmapped = pair(name)
        for read in unmapped:
            read.is_unmapped = read.mate_is_unmapped = True
            read.reference_id = read.next_reference_id = -1
            read.reference_start = read.next_reference_start = -1
            read.cigarstring = None
        records += unmapped
    single = pair("single")[0]
    single.is_paired = single.is_read1 = False
    records.append(single)
    output, _ = run(tmp_path, records, kind)
    counts = Counter(r.query_name for r in output)
    assert counts == (
        {"mapped": 2}
        if kind == "linear"
        else {"mapped": 2, "unmapped1": 2, "unmapped2": 2, "single": 1}
    )


def test_legacy_distinct_second_endpoint_preserves_both_complete_pairs(tmp_path):
    # Equal R1 starts must not discard R1 independently of the distinct R2s.
    first, second = pair("first"), pair("second", end=350)
    output, _ = run(tmp_path, first + second, "legacy")
    assert Counter(r.query_name for r in output) == {"first": 2, "second": 2}


def test_legacy_conflicting_name_only_barcode_map_is_rejected(tmp_path):
    path, csv = write_input(tmp_path, pair("shared"), "legacy")
    csv.write_text("shared,A\nshared,B\n")
    with pytest.raises(ValueError, match="Conflicting barcode"):
        TagDedup(str(path), str(path) + ".bai", str(csv)).tag_dedup_reads(
            True, str(tmp_path), "out"
        )
    assert not list(tmp_path.glob("out*.bam"))


def test_validator_spills_outstanding_mates_and_detects_reused_completed_identity():
    first, second = pair("first"), pair("second")
    with PrimaryPairValidator(max_pending=1) as validator:
        scratch = validator.index_path
        for read in [first[0], second[0], first[1], second[1]]:
            validator.observe(read, template_identity(read))
            assert len(validator.pending) <= 1
        assert validator.finish() == 2
        validator.observe(first[0], template_identity(first[0]))
        with pytest.raises(ValueError, match="repeated primary pair"):
            validator.observe(first[1], template_identity(first[1]))
    assert not scratch.exists()


def test_validator_counts_missing_mates_in_memory_and_on_disk():
    with PrimaryPairValidator(max_pending=1) as validator:
        scratch = validator.index_path
        for read in [pair("first")[0], pair("second")[0]]:
            validator.observe(read, template_identity(read))
        assert validator.spilled_count == 1
        with pytest.raises(ValueError, match="2 template identities lack a mate"):
            validator.finish()
    assert not scratch.exists()


@pytest.mark.parametrize("kind", ["linear", "legacy"])
def test_discordant_cross_contig_pair_is_kept_whole(tmp_path, kind):
    records = pair("cross_contig")
    records[0].next_reference_id = 1
    records[1].reference_id = 1
    for read in records:
        read.template_length = 0
    output, _ = run(tmp_path, records, kind)
    assert len(output) == 2
    assert {(read.is_read1, read.reference_id) for read in output} == {(True, 0), (False, 1)}


@pytest.mark.parametrize("kind", ["linear", "legacy"])
def test_partly_unmapped_pair_may_carry_mapped_mate_coordinates(tmp_path, kind):
    records = pair("partial")
    records[0].mate_is_unmapped = True
    records[0].next_reference_start = records[0].reference_start
    records[1].is_unmapped = True
    records[1].reference_start = records[0].reference_start
    records[1].cigarstring = None
    for read in records:
        read.template_length = 0
    output, _ = run(tmp_path, records, kind)
    assert len(output) == (0 if kind == "linear" else 2)
