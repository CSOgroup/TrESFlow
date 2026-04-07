# TrESFlow: Usage

## Overview

`TrESFlow` supports one public input contract: a single hierarchical YAML samplesheet passed with `--samplesheet`.
There is no CSV input mode in this repository.

The pipeline runs two independent modality branches from that YAML:

- `rna`: sample-barcode tagging, UMI tagging, cell-barcode tagging, trimming, split by SB groups, `FqToSAM`, STARsolo, filtered BAM, bigWigs
- `dna`: sample-barcode tagging, modality tagging, cell-barcode tagging, trimming, split by SB groups and DNA marks, alignment, duplicate marking, NoDup BAM, bigWig

## Quick Start

Smoke test with the bundled example YAML:

```bash
cd /mnt/pdata/nikita/TrESFlow
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Canonical real-data style run:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.real.example.yaml \
  --outdir results/test_real
```

The pipeline expects the runtime toolchain under `--runtime_env_prefix` and validates it before starting the workflow.

## Samplesheet Contract

The supported YAML structure is:

```yaml
library_name: MY_LIBRARY

resources:
  ligation_barcode_whitelist: /path/to/ligation_barcode_whitelist.txt
  rna_ref_base_dir: /path/to/reference_base_dir
  rna_align_species: human
  dna_bwa_reference: /path/to/bwa_index_prefix
  dna_blacklist_bed: /path/to/blacklist.bed
  dna_effective_genome_size: 2913022398

samples:
  day15:
    groups:
      Normal:
        sb_barcodes: [CAGT, ACGT]
      Co2:
        sb_barcodes: [GTCA, TGCA]

    rna:
      reads:
        i1: /path/to/day15_RNA_I1.fastq.gz
        r1: /path/to/day15_RNA_R1.fastq.gz
        r2: /path/to/day15_RNA_R2.fastq.gz

    dna:
      reads:
        i1: /path/to/day15_DNA_I1.fastq.gz
        i2: /path/to/day15_DNA_I2.fastq.gz
        r1: /path/to/day15_DNA_R1.fastq.gz
        r2: /path/to/day15_DNA_R2.fastq.gz
      mark_barcodes:
        H3K27me3: AGGCTATA
        H3K27ac: GCCTCTAT
```

### Top-level fields

- `library_name`: run-level library label propagated into RG headers and derived contract files
- `resources`: shared runtime resources for the whole run
- `samples`: biological sample blocks keyed by user-defined sample ID

### `resources`

These values can be defined in the YAML and, when needed, overridden by CLI params.

- `ligation_barcode_whitelist`: ligation barcode whitelist used by RNA and DNA cell-barcode tagging
- `rna_ref_base_dir`: STARsolo reference base directory
- `rna_align_species`: `human` or `mouse`
- `dna_bwa_reference`: bwa-mem2 index prefix
- `dna_blacklist_bed`: BED file for blacklist filtering
- `dna_effective_genome_size`: effective genome size used for DNA coverage

### `samples.<sample_id>.groups`

`groups` is the source of truth for biological sample-barcode grouping.

- each group key is the biological label that will appear in split outputs
- `sb_barcodes` is the list of sample barcodes assigned to that group
- sample barcodes must be unique within a sample block

### `samples.<sample_id>.rna`

The RNA block is optional, but if present it must contain:

- `reads.i1`
- `reads.r1`
- `reads.r2`

### `samples.<sample_id>.dna`

The DNA block is optional, but if present it must contain:

- `reads.i1`
- `reads.i2`
- `reads.r1`
- `reads.r2`
- `mark_barcodes`

`mark_barcodes` maps biological mark labels, such as `H3K27ac`, to their DNA modality barcodes.

## Derived Internal Contract

The parser in [`lib/SamplesheetParser.groovy`](../lib/SamplesheetParser.groovy) turns the hierarchical YAML into modality-specific work rows and writes helper contract files under:

```text
<outdir>/pipeline_info/derived_contract/
```

These derived files include:

- `sb_group_map.tsv`
- `dna_mo_map.tsv` when DNA is present
- per-sample DNA modality whitelist files

This keeps the public input contract user-friendly while preserving the split and alignment interfaces used by the current core modules.

## Parameters

The main public parameters are:

- `--samplesheet`
- `--outdir`
- `--runtime_env_prefix`
- `--max_cpus`

Optional resource override parameters remain available:

- `--ligation_barcode_whitelist`
- `--rna_ref_base_dir`
- `--rna_align_species`
- `--dna_bwa_reference`
- `--dna_blacklist_bed`
- `--dna_effective_genome_size`

## Bundled Examples

- smoke-test YAML: [`assets/samplesheet.example.yaml`](../assets/samplesheet.example.yaml)
- real-data style example: [`assets/samplesheet.real.example.yaml`](../assets/samplesheet.real.example.yaml)
- editable template: [`assets/samplesheet.template.yaml`](../assets/samplesheet.template.yaml)
