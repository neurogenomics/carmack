# Barcode edit-distance validation, October 2026

The exact bit-vector implementation at `640fe9ea1ef00eb7f8b0d05125b3d6b48a78e98a`
was compared with its parent, `9dc74d89fb26e3bd11459d8183836631a2417bdc`.
It preserves global unit-cost Levenshtein distance, including symmetric `N`
wildcards. Candidate windows, chemistry, allowed errors and ambiguity ranking
are unchanged.

## Correctness checks

The candidate passed the hosted container suite with 2,135 tests and 58 skips,
and all 87 golden-output cases passed without updating expected outputs. An
independent scalar oracle checked 266,152 distance cases and 3,000 complete
best-window tuples, including long strings, Unicode and custom wildcard markers.
The hosted result is recorded in [this Actions run](https://github.com/neurogenomics/carmack/actions/runs/38084731403).

The cluster experiment then compared every scientific product from 40 extraction
commands: eight 200-read runtime gates, eight 20,000-read warmups and 24 measured
20,000-read invocations. All corresponding baseline/candidate outputs matched.
The checks covered ordered decompressed annotated FASTQ and barcode tables,
barcode counts, MultiQC JSON, and the PNG plot. Exactly one validated generation
timestamp line was removed from the text-report comparison; its original bytes
were retained and hashed. Gzip products were compared after decompression.
Input IDs, mate synchronization, accepted-read order, sequence/quality preservation
and barcode-count reconciliation were checked separately.

This establishes computational equivalence on the tested inputs. It does not
establish biological barcode-assignment accuracy.

## Controlled cluster measurements

One CX3 node reserved 32 CPUs and 12 GB. Every command ran sequentially in the
same frozen container: Python 3.12.14, pigz 2.8, gzip 1.13 and identical libraries.
The runtime verified its 32-CPU affinity mask and the exact source overlay.
Worker counts below refer to extraction workers; three output compressors each
use four pigz threads, so they are not total CPU allocation or billing figures.

Two private paired-read workloads were selected deterministically:

- **SK609, nonprimd chemistry:** 10,000 pairs from the first 100,000 synchronized
  pairs of each of two lanes, selected by read-ID hash and restored to original
  within-lane order. This is a lane-balanced prefix workload.
- **SK508, primd chemistry:** 20,000 pairs selected by read-ID hash from an existing
  100,000-pair failure-stratified subset, retaining input order. This deliberately
  overrepresents some historical failure categories.

Only R1 is processed by barcode extraction. R2 was retained and checked for exact
pair identity. Neither workload is an unweighted representative sample of its
whole library.

After the gates and warmups, each sample/worker group had three paired measured
repetitions, with deterministic randomized group order and alternating source
versions. Wall time includes container/interpreter startup, imports, extraction,
compression and plotting. Subset preparation, full input hashing and post-run
validation are outside the per-command timing.

| Workload | Workers | Parent median (s) | Candidate median (s) | Median paired speedup | Parent / candidate wall CV |
|---|---:|---:|---:|---:|---:|
| SK609 prefix | 8 | 50.230 | 13.685 | 3.659× | 0.26% / 1.12% |
| SK609 prefix | 16 | 27.602 | 8.629 | 3.199× | 0.52% / 0.66% |
| SK508 failure-stratified | 8 | 225.590 | 57.537 | 3.894× | 1.90% / 0.49% |
| SK508 failure-stratified | 16 | 116.108 | 30.805 | 3.769× | 0.69% / 0.47% |

Speedup is the median of the three within-pair wall-time ratios. It need not equal
the ratio of the two wall-time medians. These are warm, small-workload results
from one node. They do not establish full-library or full-pipeline speedup,
cohort-wide throughput, scheduler efficiency or an optimal worker count.

Source, container, complete input and selected-subset hashes were pinned. Every
run retained raw product hashes, scientific comparison hashes, runtime identity
and GNU-time CPU/RSS evidence. The final publication check rehashed all products
and closed the source/output identity windows. An independent collection review
verified 411 small artifacts and recalculated all 12 paired ratios. Large product
files stayed on the cluster; their equivalence is supported by the reviewed
runtime and its complete final binding, rather than a second local full-file scan.

Evidence identities:

- Runtime result SHA-256: `b8ba7c5a9e9d8fddc3a2b7e4dacb0553fcc9072e93fc5f5c38ed34fc01f0f3d1`.
- Collection receipt SHA-256: `d1f7baead9da5b5caf6e4d78f0e62bb10c1fbbff1b0c14d1e571eb3137e2a0ef`.
- Reviewed execution-plan SHA-256: `3369342322e15ba4685a4c184f8504d4eb2fbafa2dbe857f9844c3767f2c2e7e`.
- Container SHA-256: `ccc4ce0e88addee2cee1a038c279f628e3f832ee226d7e9ac16e03f7ec2167b7`.

The original 2,000-read local SK661 fixture also produced identical pinned
outputs: median wall time changed from 35.630 to 17.278 seconds with one worker,
and 14.750 to 8.309 seconds with four. Those local runs used gzip, ran parent
repetitions before candidate repetitions and had substantial shared-machine
variation. Keep them separate from the controlled cluster measurements above.
