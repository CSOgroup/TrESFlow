# TrESFlow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - 2026-08-24

### Added

- Multi-FASTQ input support for ordered technical sequencing chunks using YAML sequences or comma-separated paths.
- Input FASTQ provenance tracking for multi-file libraries.

### Changed

- DNA read groups distinguish AVITI run, flowcell, and lane while preserving the logical library across sequencing units.

### Fixed

- Multi-FASTQ samplesheet normalization and validation for YAML lists and comma-separated inputs.

## v1.0.0 - 2026-08-18

First public TrESFlow release.

### Added

- Nextflow DSL2 workflow for joint TrES-seq RNA and DNA preprocessing.
- Hierarchical YAML samplesheet supporting RNA-only, DNA-only, and multimodal samples.
- Single- and dual-tagmentation DNA library support.
- RNA STARsolo alignment and filtered-cell BAM/coverage outputs.
- DNA alignment, duplicate marking, NoDup BAMs, and BigWig coverage tracks.
- AVITI-aware optical duplicate handling with lane-level read groups.
- Dual-tagmentation residual-linker artifact filtering.
- FastQC, samtools QC, MultiQC, and the self-contained TrESFlow HTML QC report.
- Read-retention, barcode-composition, library-complexity, and sequencing metrics.
- Configurable resource profiles and runtime/work-directory handling.

### Changed

- Split FASTQ publication is optional while internal uncompressed FASTQs are used downstream.
- RNA filtered BAM processing and canonical BAM validation were optimized to reduce redundant I/O.
- RNA coverage generation was parallelized without changing coverage definitions.
- Output structure was consolidated under `rna_align/`, `dna_align/`, `TrES_Stats/`, and `pipeline_info/`.

### Fixed

- RNA and DNA split-retention metric propagation.
- RNA group routing for multiple sample-barcode groups.
- DNA ligation-index handling for single versus dual tagmentation.
- Picard MarkDuplicates performance regression caused by excessive read-group cardinality.
- Canonical chromosome filtering and NoDup coverage consistency.
