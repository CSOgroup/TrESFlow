# FASTQ preprocessing fixture

These two deterministic FASTQ files are synthetic, test-only inputs for the
Phase 3 Trim Galore, Cutadapt, and pigz process tests. They were authored for
TrESFlow and contain no external biological source sequence.

The `salvage` pair contains an Illumina adapter followed by an audited TrES
dual-tag signature, so ordinary trimming retains its genomic prefix. The
`residual` pair contains an internal exact signature and is rejected by the
dual-tag filter. The `clean` pair is retained unchanged. The files are covered
by the repository MIT license and their SHA-256 hashes are recorded in
`manifest.json`.
