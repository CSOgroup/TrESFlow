# Portable Python helper fixtures

These files are deterministic, synthetic, test-only inputs for the isolated
`BARCODE_GATE_METRICS` and `TRES_REPORT_HTML` process tests. They contain no
human or other source biological sequence: the five paired FASTQ records use
repeated artificial `ACGT`/`TGCA` sequences, and the TSV values were authored
to exercise the nested barcode gates and two report branches.

The fixtures were authored for TrESFlow Phase 1, are distributed under the
repository's MIT license, and have no external source or download dependency.
Their byte-level provenance is recorded in `manifest.json`; the Python test
suite verifies that every listed SHA-256 digest matches the committed file.

Expected gate progression is 5 input pairs, 4 ligation-accepted pairs, 3
sample-barcode-accepted pairs, and 2 modality-barcode-accepted pairs. The
accepted branches are `run-alpha/mark-one` and `run-beta/mark-three`.
