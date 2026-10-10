# Cell-calling input and output contract

```sh
mkdir -p cells
MPLBACKEND=Agg carmack call-cells peaks.bed aligned.tagged.bam \
  --force_n 2000 --min_overlap 1 --visualise --prefix sample --output_dir cells
```

The BAM must be coordinate sorted, indexed and carry a `CB` tag on every alignment read
by the command. Supply an explicit BAI path as the third positional argument when it is
not named `<bam>.bai`. Contig names must match the BED. The current BED reader requires
at least six columns; use unique names in column 4 and an integer score in 0–1000 in
column 5. Bare BED3 and comment/header lines are not supported.

Each matrix entry counts **alignment records** with at least `min_overlap` aligned query
bases inside a peak. CIGAR deletions and reference skips do not supply aligned bases.
Overlapping mates can contribute two counts to the same peak, and one alignment can
contribute to multiple overlapping peaks. The caller does not itself remove duplicates,
secondary/supplementary alignments or low-MAPQ records. Apply and record those filters
upstream. These are not automatically unique-fragment, molecule or insertion-site counts.

Without `--force_n`, a decreasing concave knee heuristic chooses how many barcode columns
to retain. This is a heuristic needing experimental validation. A flat/singleton/empty
curve or an absent knee now raises an actionable error instead of silently selecting
every barcode. `--force_n 0` explicitly retains all barcodes, including an empty BAM;
an empty rank curve cannot be plotted. Tied counts are ordered by barcode alphabetically.

Output files are `sample_matrix.mtx`, `sample_peaks.bed`, `sample_barcodes.tsv`, and,
when requested, `sample_barcode_matrix.png`. Matrix rows correspond **in order** to
exported peak coordinates; columns correspond in order to exported barcodes. Only peaks
with nonzero counts among selected barcodes are emitted. Exact repeated genomic intervals
are counted and exported once, using the first interval record. Different intervals
cannot share a peak name. Peak rows are sorted by name, so the exported BED need not be
coordinate sorted; sorting that BED alone would break its correspondence with the matrix.

Validate exports by checking both axis lengths and named counts. Matching dimensions
alone will not reveal a swapped row label. The review regression fixtures use known
counts at deliberately nonalphabetical peak coordinates and explicit CIGAR gaps. These
establish software behavior; they do not validate cell identities on real multiome data.
