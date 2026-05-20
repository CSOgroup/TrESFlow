# TrESFlow: Usage

## Overview

`TrESFlow` supports one public input contract: a single hierarchical YAML samplesheet passed with `--samplesheet`.
There is no CSV input mode in this repository.

The pipeline runs two independent modality branches from that YAML, then builds sequencing-efficiency UpSet plots from the published tag-record and alignment channels:

- `rna`: sample-barcode tagging, UMI tagging, cell-barcode tagging, trimming, split by SB groups, `FqToSAM`, STARsolo, filtered BAM, bigWigs
- `dna`: sample-barcode tagging, modality tagging, cell-barcode tagging, trimming, split by SB groups and DNA marks, alignment, duplicate marking, NoDup BAM, bigWig

Sequencing-efficiency outputs are written to `TrES_Stats/` as UpSet PDFs only. Sankey plots, HTML reports, count tables, combined RNA+DNA reports, and sequencing-efficiency warning TSVs are not produced. Optional unavailable BAM-derived categories are skipped with warnings in the process log.

DNA alignment no longer filters out low-count cell barcodes during `ALIGN_DNA`. Low-count status is visualized in sequencing-efficiency plots as `CB>100 +`, using BAM-derived unique query-name read-pair support.

## Quick Start

Smoke test with the bundled example YAML:

```bash
cd /mnt/dataFast/ahrmad/tresflowdir/TrESFlow
NXF_OFFLINE=true nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Canonical real-data style run:

```bash
NXF_OFFLINE=true nextflow run . \
  --samplesheet /mnt/dataFast/ahrmad/TEST_NF/isa_multiome.yaml \
  --outdir /mnt/dataFast/ahrmad/TEST_NF/TrESFlow_Isa \
  --max_cpus 32
```

The pipeline reads runtime and reference locations from the samplesheet. Runtime and reference CLI overrides are rejected.

## Samplesheet Contract

The supported YAML structure is:

```yaml
library_name: Isa

runtime:
  env_prefix: /home/annan/micromamba/envs/tres

references:
  species: human
  root: /mnt/dataFast/ahrmad/TrESFlow_References
  ligation_barcode_whitelist: /mnt/dataFast/ahrmad/TrESFlow_References/ligation_barcode_whitelist.txt
  rna_ref_dir: /mnt/dataFast/ahrmad/TrESFlow_References/rna/human/star
  dna_ref_dir: /mnt/dataFast/ahrmad/TrESFlow_References/dna/human/bwa
  dna_blacklist_bed: /mnt/dataFast/ahrmad/TrESFlow_References/dna/human/hg38-blacklist.v2.bed
  dna_chrom_sizes: /mnt/dataFast/ahrmad/TrESFlow_References/dna/human/hg38.chrom.sizes
  dna_effective_genome_size: 2913022398

samples:
  day15:
    groups:
      Normal:
        rna_sb_barcodes: [CAGT, ACGT]
        dna_sb_barcodes: [CAG, ACG]
      Co2:
        rna_sb_barcodes: [GTCA, TGCA]
        dna_sb_barcodes: [GTC, TGC]

    rna:
      reads:
        i1: /path/to/day15_RNA_I1.fastq.gz
        r1: /path/to/day15_RNA_R1.fastq.gz
        r2: /path/to/day15_RNA_R2.fastq.gz

    dna:
      tagmentation: dual
      reads:
        i1: /path/to/day15_DNA_I1.fastq.gz
        r1: /path/to/day15_DNA_R1.fastq.gz
        r2: /path/to/day15_DNA_R2.fastq.gz
      mark_barcodes:
        H3K27me3: AGGCTATA
        H3K27ac: GCCTCTAT
```

### Top-level fields

- `library_name`: run-level library label propagated into RG headers and derived contract files
- `runtime`: required runtime environment and explicit task temporary directory
- `references`: required species label, shared files, and direct RNA/DNA reference paths
- `samples`: biological sample blocks keyed by user-defined sample ID

### `runtime`

- `env_prefix`: environment prefix containing `python3`, `codon`, `trim_galore`, `STAR`, `samtools`, `bedGraphToBigWig`, `bwa-mem2`, `bamCoverage`, and `gatk`
- `tmpdir`: optional explicit task temporary directory. If omitted, the pipeline uses `--outdir`. The pipeline creates it if missing and fails if it is not writable.

### `references`

`references.rna_ref_dir` points directly to the STAR index directory. The pipeline passes this exact path to STAR and does not append species, `rna`, or `star`.

`references.dna_ref_dir` points to the directory containing exactly one complete bwa-mem2 sidecar set. The inferred prefix is used for `bwa-mem2 mem`.

`references.dna_effective_genome_size` is required for DNA runs because `BAM_COVERAGE_DNA` passes it to `bamCoverage --effectiveGenomeSize`.

### `samples.<sample_id>.groups`

`groups` is the source of truth for biological sample-barcode grouping.

- each group key is the biological label that will appear in split outputs
- `rna_sb_barcodes` and `dna_sb_barcodes` are modality-specific sample barcodes assigned to that logical group.
- Legacy `sb_barcodes` remains supported for single-tagmentation samples. In `dna.tagmentation: dual`, DNA requires explicit 3 nt `dna_sb_barcodes`; they are not derived from RNA barcodes.
- sample barcodes must be unique within a sample block

### `samples.<sample_id>.rna`

The RNA block is optional, but if present it must contain:

- `reads.i1`
- `reads.r1`
- `reads.r2`

### `samples.<sample_id>.dna`

The DNA block is optional, but if present it must contain:

- `reads.i1`
- `reads.r1`
- `reads.r2`
- `mark_barcodes`

`reads.i2` is required for `dna.tagmentation: single` and optional for `dna.tagmentation: dual`.

DNA ligation tagging uses the same `Tag_Lig3` correction and output format for both modes, but with mode-specific barcode-source reads and start positions:

- `single`: ligation source `reads.i1`, L1/L2/L3 starts `15,53,91`
- `dual`: ligation source `reads.i1`, L1/L2/L3 starts `41,79,117`

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
- `--max_cpus`
- `--cleanup_work`
- `--rna_starsolo_cpus`
- `--dna_align_cpus`
- `--coverage_cpus`
- `--helper_cpus`
- `--tagging_cpus`
- `--tagging_memory`

Deprecated runtime/reference CLI parameters now fail with a hard error.

For local execution, `--max_cpus` is the global executor cap and all bundled per-process CPU reservations are capped by it. The default reservations favor concurrency across independent samples, groups, and DNA marks:

- RNA STARsolo and DNA alignment default to `16` CPUs each.
- RNA and DNA coverage default to `8` CPUs.
- trim, split, RNA filtered-BAM, and DNA duplicate-filter helpers default to `4` CPUs.
- barcode-tagging steps default to `4` CPUs and `32 GB` memory.
- `FQ_TO_SAM` and `MARK_DUPLICATES_DNA` stay at `1` CPU.

Override the bucket params above on the command line or in a Nextflow config when a specific machine or scheduler profile can support larger reservations.

`--cleanup_work` defaults to `true`. TrESFlow uses Nextflow's supported successful-run cleanup to remove task work directories after final outputs are published and all downstream consumers have completed. This keeps large FASTQ, uSAM, tag-record, and BAM intermediates from remaining in `work/` after a successful run. The tradeoff is that `--resume` is not expected to be reliable for cleaned tasks. Set `--cleanup_work false` for debugging or for runs where preserving work directories is more important than disk cleanup.

## Bundled Examples

- smoke-test YAML: [`assets/samplesheet.example.yaml`](../assets/samplesheet.example.yaml)
- real-data style example: [`assets/samplesheet.real.example.yaml`](../assets/samplesheet.real.example.yaml)
- editable template: [`assets/samplesheet.template.yaml`](../assets/samplesheet.template.yaml)
