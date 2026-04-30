# TrESFlow: Output

## Overview

`TrESFlow` writes modality-specific outputs directly under `--outdir`.
The exact directories present depend on whether the YAML samplesheet contains `rna`, `dna`, or both modalities.

## Top-level result directories

RNA-related outputs:

- `split/`
- `align/`

DNA-related outputs:

- `dna_split/`
- `dna_align/`

Shared reporting outputs:

- `pipeline_info/`

## RNA outputs

### `split/`

Per-group RNA FASTQs and RG headers from the split stage:

- `<sample>_<group>_R1.fq.gz`
- `<sample>_<group>_R2.fq.gz`
- `SAM_RG_Header_<sample>_<group>.tsv`

### `align/`

STARsolo and filtered BAM outputs:

- `<sample>_<group>.Solo.outGeneFull/`
- `<sample>_<group>.filtered_cells.bam`
- `<sample>_<group>.stranded_*.bw`
- `<sample>_<group>.unstranded_*.bw`

## DNA outputs

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

- published RNA FASTQs are the grouped split FASTQs under `split/`
- published DNA FASTQs are the grouped and marked split FASTQs under `dna_split/`
- earlier tag, trim, RNA uSAM, STAR aligned BAM, duplicate-marking, NoDup BAM, and non-published coverage files remain transient in `work/` unless captured manually
