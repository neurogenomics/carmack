# Read processing walkthrough

Carmack prepares reads and provides BAM utilities. Genome alignment, reference preparation,
UMI correction and a complete multiome analysis are separate steps. Record the Carmack Git
commit or container digest, chemistry, input checksums, reference build, command lines and
dependency versions with each analysis.

## Choose the chemistry

| Name | R1 layout, 5′ to 3′ | Supported read stages |
|---|---|---|
| `carmack_custom_seq_1_0` | BC3(10), primer C(22), BC2(10), primer A(22), BC1(10), UMI(8), poly-G, target index(8), mosaic end(19) | All four below |
| `carmack_custom_seq_1_0_primd` | Leading primer D(22), then the same layout | All four below |
| `hydrop` | BC3(10), spacer(10), BC2(10), spacer(10), BC1(10) | Barcode extraction |

Lengths are bases. The leading primer D's sequence is unspecified. Barcode order is
**BC3 + BC2 + BC1**, matching read order, in both output arms. The shipped barcode and
target whitelists live under `carmack/data/`. A matching chemistry name is not evidence
that a library has that layout: verify the sample sheet and inspect positional and
homopolymer diagnostics on a small representative subset first.

## Run the four read stages

Create the output directory before running commands. Use an explicit, unique prefix so
dot-containing input names and repeated runs cannot accidentally share an output name.
For batch execution and on macOS, `MPLBACKEND=Agg` avoids an interactive plotting backend.

```sh
export MPLBACKEND=Agg
mkdir -p results
carmack extract-barcodes -c carmack_custom_seq_1_0 -n 4 -p sample -o results reads_R1.fastq.gz
carmack extract-umis -c carmack_custom_seq_1_0 -p sample -o results results/sample.r1_annotated.fastq.gz
carmack assign-targets -c carmack_custom_seq_1_0 -n 4 -p sample -o results results/sample.r1_umi.fastq.gz
carmack prepare-reads -c carmack_custom_seq_1_0 -n 4 -p sample -o results results/sample.r1_tgidx.fastq.gz reads_R2.fastq.gz
```

R1 and R2 must derive from the same paired library. R1 survivors remain in input order;
`prepare-reads` searches forward through the original R2 stream by the first whitespace
token of the read name. That token must match **exactly** between mates: `/1` and `/2`
suffixes are not normalized. Reordering or independently sampling the mates breaks this
contract. Use consistent unique read names when combining lanes or libraries.

The parser accepts four-line plain, gzip and LZ4 FASTQ (LZ4 needs the `lz4` executable).
It checks non-empty `@` headers, `+` separators, matching sequence/quality lengths and
complete records. It does not support wrapped/multiline FASTQ. Framing checks do not
establish sample identity or correct library chemistry. A failed stage may leave partial
output files: accept outputs only after a zero exit status and successful reconciliation,
and rerun in a clean directory after failure.

## Check the outputs

| Stage | Main handoff | What to check |
|---|---|---|
| `extract-barcodes` | `sample.r1_annotated.fastq.gz` | Accepted + rejected = raw R1 count; barcode rank, correction and ambiguity rates |
| `extract-umis` | `sample.r1_umi.fastq.gz` | Accepted + missing-anchor + truncated = input; UMI width and homopolymer distribution |
| `assign-targets` | `sample.r1_tgidx.fastq.gz` | Output count = input; matched + no-anchor + short-window + no-match = input |
| `prepare-reads` | Per-target R1/R2 and unmatched R1/R2/barcodes triples | Written + insert-not-sequenced = input; equal counts and IDs within each output pair/triple |

Each stage emits a text stats file and MultiQC custom-content JSON files. See
[prepare-reads outputs](prepare-reads-outputs.md) for exact final filenames. Use
`sample.detected_targets.txt` to find populated buckets; empty gzip files can exist for
unused targets. Read `sample.none.barcodes.json` for barcode/UMI positions instead of
hardcoding them in downstream tools. Positions in annotation tags and the sidecar are
0-based half-open spans.

Barcode `PERFECT`, `CORROK` and `FAIL` are mutually exclusive read outcomes. Component
ambiguity is a subset of failures, not a fourth disjoint class. The MultiQC field named
`pct_ambiguous` is retained for compatibility but is labelled **ambiguity events per
100 reads**: summing BC1/BC2/BC3 ambiguity can exceed 100 when a read has several ambiguous
components. Derive a distinct ambiguous-read fraction from `bc_all.txt.gz` by counting
each read with any component `-AMBIG` once. Target `NONE` currently combines distinct
unassigned reasons; use `tgidx_stats.txt` for the mutually exclusive aggregate counts.

The scTIP QNAME is `read_id|CB=<barcode>|UR=<umi>`; an aligner preserving QNAME does **not**
automatically turn those strings into BAM `CB`/`UR` tags. A downstream conversion step
must do that before tools requiring tags. The scRNA arm emits a corrected barcode plus
the raw UMI in a separate barcodes FASTQ; consult the aligner's own input contract when
using the R1/R2/barcodes triple.

## Scientific interpretation

UMIs are fixed-width raw slices; extraction does not correct sequencing errors or count
molecules. `TGIDX=NONE` means no target was assigned. It can include true scRNA molecules
and damaged/unrecognized scTIP reads; it is a routing decision, not ground-truth modality.
The barcode whitelist has a known BC2 pair, `AGCTTGAGAG` and `GGCTTGAGAG`, separated by one
substitution. An error turning one into the other is indistinguishable from a perfect
barcode. HyDrop has additional close pairs under its larger correction budget. Changing
those whitelists requires a protocol/version decision for existing libraries.

`linear-dedup` groups paired primary R1 alignments by cell barcode, chromosome, strand and
R1 start (forward) or end (reverse), keeping the highest R1 `AS` score. It is a **one-end
position rule**, not UMI deduplication or a two-end fragment identity rule. Reads without
`AS` have lowest priority; exact ties keep the first R1 encountered. Both primary mates
must carry the configured cell-barcode tag. Template selection uses **(QNAME, RG, cell
barcode)**, so a name reused in another cell or read group cannot resurrect a discarded
pair. An absent RG is its own namespace. RG distinguishes template names; it does not
split the biological duplicate groups, which still compare lanes within a cell.

Before output is published, every eligible mapped primary template must have exactly
one R1 and one R2 with reciprocal mate coordinates and strand flags. Missing mates,
inconsistent cell/RG tags, and multiple primary records in the same identity namespace
raise an error naming an example and, for incomplete pairs, the number of affected
identities. Carmack cannot reconstruct physical pairing for genuinely indistinguishable
reused identifiers and does not guess from record order. Secondary/supplementary,
unpaired and unmapped-pair records remain excluded from `linear-dedup`. The final output
reconciles to two records per winning template. Invalid input does not publish a new BAM
or successful stats report. Validate the one-end collision model against the assay
before interpreting output as molecule counts.

The older `bam-tag-deduplicate` command consumes a two-column read-ID/barcode CSV, which
is not the current annotated `bc_valid.txt.gz` format. Its optional UMI map is headerless
TSV: read ID, barcode, raw UMI, corrected UMI. It now validates complete primary pairs
and elects the first primary R1 for each **cell + ordered pair of reference starts,
strands and absolute template length**, adding corrected UMI when a map is supplied.
The resulting DU decision applies to both mates and their secondary/supplementary
alignments. Sharing only one endpoint is insufficient to discard either mate. This is
deliberately different from the linear command's one-end/AS rule.

Legacy unpaired records and unmapped pairs are retained without positional collapse;
fully unmapped records are included in both BAM and record-count accounting. A
secondary/supplementary record with no primary decision is retained uncollapsed. Records
absent from the barcode CSV are still skipped and logged. The name-only CSV cannot
express different barcodes for the same QNAME: conflicting barcode/UMI-map assignments
are rejected rather than silently overwritten. Missing UMI-map entries retain the
previous fallback to position-only grouping among entries without a corrected UMI.
Legacy duplicate stats count alignment records, including unmapped records, not molecules.

Pair validation uses a temporary on-disk SQLite uniqueness index, with an 8 MiB page-cache
target and at most 100,000 outstanding mates in memory; older unmatched mates spill to
the same database. Set `TMPDIR` to writable job scratch with sufficient capacity. Its
disk use grows with eligible template count and is logged at successful validation;
it is removed on normal completion or a handled error. An abrupt process kill can leave
scratch behind. This avoids an additional whole-BAM identity set in RAM. The linear winner set still scales with kept templates and its duplicate-group
accumulator with the largest chromosome. Benchmark peak RSS and scratch use on
representative full-depth BAMs before changing job memory requests.
