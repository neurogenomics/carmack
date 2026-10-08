# prepare-reads outputs

`prepare-reads` writes everything into `--output_dir` (default `.`). Every file name starts
with `<prefix>`, which is `--prefix` if given and otherwise the R1 input's file name up to
its first `.`.

| File | Contents |
|---|---|
| `<prefix>.none.r1.fastq.gz` | scRNA arm R1: reads with no target index, trimmed to the cDNA insert. |
| `<prefix>.none.r2.fastq.gz` | scRNA arm R2, passed through unchanged. |
| `<prefix>.none.barcodes.fastq.gz` | scRNA arm barcodes record: the corrected cell barcode followed by the UMI. |
| `<prefix>.none.barcodes.json` | Layout of the barcodes record (see below). |
| `<prefix>.<TGIDX>.r1.fastq.gz`, `<prefix>.<TGIDX>.r2.fastq.gz` | One scTIP pair per entry in the chemistry's target index whitelist. R1 is trimmed past the target index; both reads carry the same pipe-delimited name with the barcode and UMI. |
| `<prefix>.prepare_stats.txt` | Human-readable, reconciling report of the run. |
| `<prefix>.detected_targets.txt` | `token<TAB>count` lines, no header, for each arm that received at least one read: `NONE` for the scRNA arm, otherwise the target index. |
| `<prefix>.prepare_general_stats_mqc.json` | MultiQC general stats payload. |
| `<prefix>.prepare_target_distribution_mqc.json` | MultiQC per-target read distribution payload. |

Every whitelisted target gets its pair of files, even if no read lands in it, so use
`detected_targets.txt` rather than globbing the directory to find the targets that have
reads. A chemistry with no target index writes no `<TGIDX>` pairs.

The FASTQs are opened at the start of the run. The JSON layout file, the two text reports
and the MultiQC payloads are written only once every read has been processed, so a run
that fails part-way leaves none of them behind. A completed run always writes them, even
for empty input.

## Barcodes layout

`<prefix>.none.barcodes.json` describes the records in `<prefix>.none.barcodes.fastq.gz`,
so a consumer can find the cell barcode and UMI without knowing the chemistry. Each record
is every `BARCODE` component's corrected value in read-structure order, then the UMI, at
their declared lengths, so every record has the same length.

| Field | Meaning |
|---|---|
| `chemistry` | Chemistry name. |
| `length` | Total length of a barcodes record. |
| `cell_barcode` | `start` and `length` of the concatenated cell barcode. |
| `umi` | `start` and `length` of the UMI. |
| `components` | Each component's `name`, `type` (`BARCODE` or `UMI`), `start` and `length`, in record order. |

Starts are 0-based offsets into the barcodes record, not into the original R1. The file is
tool-neutral: translating it into a particular tool's options, such as STARsolo's 1-based
`--soloCBstart`, is up to the consumer.

For example, for `carmack_custom_seq_1_0` (whitespace differs in the real file):

```json
{
  "chemistry": "carmack_custom_seq_1_0",
  "length": 38,
  "cell_barcode": { "start": 0, "length": 30 },
  "umi": { "start": 30, "length": 8 },
  "components": [
    { "name": "BC3", "type": "BARCODE", "start": 0,  "length": 10 },
    { "name": "BC2", "type": "BARCODE", "start": 10, "length": 10 },
    { "name": "BC1", "type": "BARCODE", "start": 20, "length": 10 },
    { "name": "UMI", "type": "UMI",     "start": 30, "length": 8 }
  ]
}
```
