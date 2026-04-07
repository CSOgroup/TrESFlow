# TrESFlow: Output

## Overview

`TrESFlow` writes modality-specific outputs directly under `--outdir`.
The exact directories present depend on whether the YAML samplesheet contains `rna`, `dna`, or both modalities.

## Top-level result directories

RNA-related outputs:

- `tagging/`
- `split/`
- `usam/`
- `align/`

DNA-related outputs:

- `dna_tagging/`
- `dna_split/`
- `dna_align/`
- `dna_dedup/`
- `dna_nodup/`
- `dna_coverage/`

Shared reporting outputs:

- `qc/`
- `pipeline_info/`

## RNA outputs

### `tagging/`

RNA tagging publishes barcode metrics and final tagged FASTQs from the RNA branch:

- `<sample>.sample_barcode.counts.tsv`
- `<sample>.sample_barcode.stats.tsv`
- `<sample>.umi.counts.tsv`
- `<sample>.cell.counts.tsv`
- `<sample>.tag_records.tsv`
- `<sample>.cell.stats_L1.tsv`
- `<sample>.cell.stats_L2.tsv`
- `<sample>.cell.stats_L3.tsv`

Intermediate tagging FASTQs are not retained beyond the final published outputs.

### `split/`

Per-group RNA FASTQs and RG headers from the split stage:

- `<sample>_<group>_R1.fq.gz`
- `<sample>_<group>_R2.fq.gz`
- `SAM_RG_Header_<sample>_<group>.tsv`

### `usam/`

Unmapped SAM files created from the grouped RNA FASTQs:

- `<sample>_<group>_tagged.usam`

### `align/`

STARsolo and filtered BAM outputs:

- `<sample>_<group>.Solo.outGeneFull/`
- `<sample>_<group>.Aligned.sortedByCoord.out.bam`
- `<sample>_<group>.filtered_cells.bam`
- `<sample>_<group>.stranded_*.bw`
- `<sample>_<group>.unstranded_*.bw`

## DNA outputs

### `dna_tagging/`

DNA tagging publishes barcode metrics and tag-record tables:

- `<sample>.dna_sample_barcode.counts.tsv`
- `<sample>.dna_sample_barcode.stats.tsv`
- `<sample>.dna_modality.counts.tsv`
- `<sample>.dna_modality.stats.tsv`
- `<sample>.dna_cell.counts.tsv`
- `<sample>.dna_tag_records.tsv`
- `<sample>.dna_cell.stats_L1.tsv`
- `<sample>.dna_cell.stats_L2.tsv`
- `<sample>.dna_cell.stats_L3.tsv`

### `dna_split/`

Per-group and per-mark DNA split outputs:

- `<sample>_<group>_<mark>_R1.fq.gz`
- `<sample>_<group>_<mark>_R2.fq.gz`
- `SAM_RG_Header_<sample>_<group>_<mark>.tsv`

### `dna_align/`

Filtered aligned BAMs and per-barcode read counts:

- `<sample>_<group>_<mark>.bam`
- `<sample>_<group>_<mark>.bam.bai`
- `<sample>_<group>_<mark>_ProperPairedMapped_reads_per_barcode.tsv`

### `dna_dedup/`

Duplicate-marked BAMs and Picard/GATK duplicate metrics:

- `<sample>_<group>_<mark>_MarkedDup.bam`
- `<sample>_<group>_<mark>_MarkedDup.bam.bai`
- `<sample>_<group>_<mark>.DuplicateMetrics.txt`

### `dna_nodup/`

Duplicate-filtered DNA BAMs:

- `<sample>_<group>_<mark>_NoDup.bam`
- `<sample>_<group>_<mark>_NoDup.bam.bai`

### `dna_coverage/`

Coverage tracks generated from the NoDup BAMs:

- `<sample>_<group>_<mark>_NoDup.bw`

## QC outputs

The `qc/` directory contains run-level summary tables and plot folders:

- `rna_sample_stage_counts.tsv`
- `rna_group_stage_counts.tsv`
- `dna_sample_stage_counts.tsv`
- `dna_group_stage_counts.tsv`
- `dna_group_mark_stage_counts.tsv`
- `qc/rna/` plots
- `qc/dna/` plots

## Pipeline information

`pipeline_info/` contains execution metadata and the derived helper contract written from the YAML samplesheet.

Expected files include:

- `execution_report.html`
- `execution_timeline.html`
- `execution_trace.tsv`
- `flowchart.html`
- `runtime_contract.tsv`

When the YAML contains group and DNA mark definitions, the parser also writes:

```text
pipeline_info/derived_contract/
```

with files such as:

- `sb_group_map.tsv`
- `dna_mo_map.tsv`
- per-sample DNA modality whitelist files

## FASTQ retention policy

The pipeline does not keep every intermediate tagging or trimming FASTQ in the published results.

- published RNA FASTQs are the grouped split FASTQs under `split/`
- published DNA FASTQs are the grouped and marked split FASTQs under `dna_split/`
- earlier tag and trim FASTQs remain transient in `work/` unless captured manually
