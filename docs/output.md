# TrESFlow: Output

## Overview

`TrESFlow` writes modality-specific outputs directly under `--outdir`.
The exact directories present depend on whether the YAML samplesheet contains `rna`, `dna`, or both modalities.

## Top-level result directories

RNA-related outputs:

- `rna_split_fastqs/`
- `rna_align/`
- `TrES_Stats/`

DNA-related outputs:

- `dna_split_fastqs/`
- `dna_align/`
- `TrES_Stats/`

Shared reporting outputs:

- `pipeline_info/`

## RNA outputs

### `rna_split_fastqs/`

Per-group RNA FASTQs from the split stage:

- `<sample>_<group>_R1.fastq.gz`
- `<sample>_<group>_R2.fastq.gz`

SAM read-group header TSVs are internal work files and are not published.

### `rna_align/`

STARsolo and filtered BAM outputs:

- `<sample>_<group>.Solo.outGeneFull/`
- `<sample>_<group>.filtered_cells.bam`
- `<sample>_<group>.stranded_*.bw`
- `<sample>_<group>.unstranded_*.bw`

## DNA outputs

### `dna_split_fastqs/`

Per-group and per-mark DNA split FASTQs:

- `<sample>_<group>_<mark>_R1.fastq.gz`
- `<sample>_<group>_<mark>_R2.fastq.gz`

SAM read-group header TSVs are internal work files and are not published.

### `dna_align/`

Final duplicate-filtered BAMs, coverage tracks, and per-barcode read counts:

- `<sample>_<group>_<mark>_NoDup.bam`
- `<sample>_<group>_<mark>_NoDup.bam.bai`
- `<sample>_<group>_<mark>_NoDup.bw`
- `<sample>_<group>_<mark>_ProperPairedMapped_reads_per_barcode.tsv`

## Tagging/count/stat outputs

### `TrES_Stats/`

Tagging summaries are published with modality-specific names:

- `<sample>.rna_sample_barcode.counts.tsv`
- `<sample>.rna_sample_barcode.stats.tsv`
- `<sample>.rna_umi.counts.tsv`
- `<sample>.rna_cell.counts.tsv`
- `<sample>.rna_cell.stats_L1.tsv`
- `<sample>.rna_cell.stats_L2.tsv`
- `<sample>.rna_cell.stats_L3.tsv`
- `<sample>.rna_tag_records.tsv.gz`
- `<sample>.dna_sample_barcode.counts.tsv`
- `<sample>.dna_sample_barcode.stats.tsv`
- `<sample>.dna_modality.counts.tsv`
- `<sample>.dna_modality.stats.tsv`
- `<sample>.dna_cell.counts.tsv`
- `<sample>.dna_cell.stats_L1.tsv`
- `<sample>.dna_cell.stats_L2.tsv`
- `<sample>.dna_cell.stats_L3.tsv`
- `<sample>.dna_tag_records.tsv.gz`

Only published tag-record tables are gzipped. Uncompressed tag-record TSVs are internal work files.

## Pipeline information

`pipeline_info/` contains execution metadata and the derived helper contract written from the YAML samplesheet.

Expected files include:

- `execution_report.html`
- `execution_timeline.html`
- `execution_trace.tsv`
- `flowchart.html`
- `runtime_contract.tsv`
- `warnings/*.zero_mapped_nodup_bam.tsv` when DNA NoDup BAMs have zero mapped reads and bamCoverage is skipped

When the YAML contains group and DNA mark definitions, the parser also writes:

```text
pipeline_info/derived_contract/
```

with files such as:

- `sb_group_map.tsv`
- `dna_mo_map.tsv`
- per-sample DNA modality whitelist files

## FASTQ retention policy

The pipeline does not keep intermediate tagging, uSAM, duplicate-marking, duplicate-split, or coverage side products in the published results.

- published RNA FASTQs are the grouped split FASTQs under `rna_split_fastqs/`
- published DNA FASTQs are the grouped and marked split FASTQs under `dna_split_fastqs/`
- earlier tag, trim, RNA uSAM, STAR aligned BAM, duplicate-marking, and non-published coverage side products remain transient in `work/` unless captured manually
