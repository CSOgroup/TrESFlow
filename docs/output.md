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

Duplicate-marked BAMs, final duplicate-filtered BAMs, and coverage tracks:

- `<sample>_<group>_<mark>_MarkedDup.bam`
- `<sample>_<group>_<mark>_MarkedDup.bam.bai`
- `<sample>_<group>_<mark>_NoDup.bam`
- `<sample>_<group>_<mark>_NoDup.bam.bai`
- `<sample>_<group>_<mark>_NoDup.bw`

`*_MarkedDup.bam` is retained so sequencing-efficiency reporting can count aligned DNA reads before duplicate removal. `*_NoDup.bam` remains the duplicate-filtered output used for downstream DNA coverage.

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

Sequencing-efficiency reports are also published under `TrES_Stats/` as UpSet PDF plots only. No Sankey plots, HTML reports, TSV/CSV count tables, combined RNA+DNA summaries, or `sequencing_efficiency.warnings.tsv` files are produced. Optional BAM-derived categories are skipped with warnings in the `SEQUENCING_EFFICIENCY` process log. Intersections are exact: the process streams category observations to task-local temporary files, sorts by unit/read id, reduces them to aggregate bitmask counts, and plots from those counts instead of storing every read id in memory.

RNA UpSet PDFs are written per sample and per sample group:

- `<sample>.rna_sequencing_efficiency.upset.pdf`
- `<sample>_<group>.rna_sequencing_efficiency.upset.pdf`

RNA UpSet categories are: `Reads +`, `Sample +`, `Ligation +`, `CB +`, `CB>100 +`, `UMI +`, `Mapped +`, and `GX +`. `Ligation +` means `L1`, `L2`, and `L3` are all present and not `NoMatch`. `CB>100 +` means the read's cell barcode has at least `--efficiency_min_read_pairs_per_cell` BAM-derived read pairs, default `100`, counted as unique query names per `CB` tag with `RG` fallback only when `CB` is absent. `GX +` requires a `GX` tag present and not `-`.

DNA UpSet PDFs are written per sample, per sample group, and per sample group plus mark:

- `<sample>.dna_sequencing_efficiency.upset.pdf`
- `<sample>_<group>.dna_sequencing_efficiency.upset.pdf`
- `<sample>_<group>_<mark>.dna_sequencing_efficiency.upset.pdf`

DNA UpSet categories are: `Reads +`, `Sample +`, `Ligation +`, `CB +`, `CB>100 +`, `Modality +`, `Mapped +`, and `Unique +`. `Mapped +` comes from mapped read names in `*_MarkedDup.bam`. `Unique +` comes from mapped read names retained in `*_NoDup.bam`, with a non-duplicate `*_MarkedDup.bam` fallback only if NoDup is unavailable or unreadable. DNA alignment no longer filters out low-count barcodes during `ALIGN_DNA`; low-count status is visualized by `CB>100 +`.

Sequencing-efficiency reports do not use `*_ProperPairedMapped_reads_per_barcode.tsv`; that file is no longer produced by `ALIGN_DNA`. These reports are not currently integrated with MultiQC.

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

The pipeline does not keep intermediate tagging, uSAM, duplicate-split, or coverage side products in the published results.

- published RNA FASTQs are the grouped split FASTQs under `rna_split_fastqs/`
- published DNA FASTQs are the grouped and marked split FASTQs under `dna_split_fastqs/`
- earlier tag, trim, RNA uSAM, STAR aligned BAM, tag-record, and non-published coverage side products are transient task outputs
- DNA duplicate-marked BAMs are published under `dna_align/` for sequencing-efficiency reporting

By default, `--cleanup_work true` asks Nextflow to clean successful task work directories after the workflow finishes successfully. This preserves published outputs while reducing retained `work/` storage. Set `--cleanup_work false` to keep work directories for debugging or a more resume-friendly run.
