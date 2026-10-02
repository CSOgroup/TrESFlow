# TrESFlow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - 2026-10-02

### Added

- Versioned physical oligo lookup (`assets/sb_oligo_lookup.v1.tsv`) with exact corresponding RNA/single-DNA and dual-DNA sample-barcode sequences, plus shared `sb_oligo_indices` selection.
- Physical indices 13–16 for RNA/single DNA: 13=`GATC`, 14=`TGAC`, 15=`CATG`, 16=`TACG`. Their dual-DNA sequences are explicitly unavailable; selecting them for dual DNA fails validation.
- Independent `rna_sb_oligo_indices` and `dna_sb_oligo_indices` selections, alongside existing sequence selectors. List order never defines RNA/DNA correspondence.
- Optional group-level `sb_oligo_pairings` to declare exceptional physical-oligo mismatches as one biological partition with a shared logical `sb_index`, while preserving each corrected physical `SB` sequence.
- Early validated identity warnings before preflight, repeated on resume, plus `pipeline_info/sb_identity_warnings.txt` and `sb_physical_to_logical.tsv` audits. Valid independent mismatches warn and continue; ambiguous or conflicting selections fail validation.

### Changed

- Dependency-aware early FASTQ cleanup under the default `--cleanup_work true`: generated FASTQs are reclaimed only after every enabled consumer succeeds, including QC and optional split-FASTQ compression. Source FASTQs and published outputs are protected; `--cleanup_work false` retains intermediates and successful task work directories.
- Complete RNA/DNA cell identity from splitting onward: `CB = XI = <sample>_<group>_<sb_index>_<L1L2L3>`, propagated through split FASTQs, BAMs, matrices, cell filtering, retention audits and per-cell statistics. Without explicit pairing, `sb_index` is the physical oligo index.
- STARsolo counts by full `XI` identity with molecular UMIs handled separately. DNA duplicate marking uses complete cell identity, preserving DNA mark routing and AVITI read-group/optical-duplicate behavior.
- Lookup-content and resolved identity changes invalidate affected resumed tasks.

### Fixed

- Prevented RNA cells with equal ligation barcodes but different sample indices from merging during counting, and prevented DNA reads from distinct cells from merging during duplicate marking.

### Compatibility and regeneration

- Existing valid samplesheet inputs remain supported, including legacy sequence selectors and modality-specific groups. Output barcode identifiers have changed from ligation-only identifiers to the complete format above.
- Previously split, aggregated or deduplicated outputs require regeneration from before splitting/counting/deduplication. Relabeling old matrices, cell calls, statistics or duplicate-marked/NoDup BAMs cannot recover merged cells or undo duplicate decisions. Use raw FASTQs or complete suitable upstream intermediates with corrected sample and ligation barcodes; normal pipeline regeneration starts from raw inputs in fresh output/work directories.

### Known limitation

- Pre-existing chromosome preflight rewrites unchanged canonical files with new modification times, so downstream filtering, coverage and reporting can rerun on an unchanged resume even when core identity tasks cache correctly.

## v1.1.1 - 2026-08-26

### Fixed

- Fixed canonical BAM filtering when noncanonical reference entries are interspersed with canonical chromosomes, preventing invalid reference IDs and BAI indexing failures.

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
