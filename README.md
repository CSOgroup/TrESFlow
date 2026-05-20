# TrESFlow

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)

TrESFlow is a Nextflow DSL2 pipeline for the preprocessing of TrES-seq data from FASTQs to cell by feature matrices.

## Install

Install your conda/mamba/micromamba env as follows (conda-forge & bioconda channels):
```bash
micromamba env create -n tres
micromamba activate tres
micromamba install pandas polars ipython pysam pybedtools numpy matplotlib seaborn scipy pyarrow upsetplot anndata scanpy matplotlib-venn leidenalg scikit-learn snapatac2
micromamba install screen samtools bwa-mem2 star fastqc multiqc trim-galore deeptools parallel ucsc-bedGraphToBigWig nextflow git gatk4
```

Download the repo and cd in it:
```bash
git clone git@github.com:CSOgroup/TrESFlow.git
cd TrESFlow
```

Install codon in your env:
```bash
./scripts/install_codon_0.16.3.sh --prefix /path/to/env/prefix
```

## Samplesheet Contract

The only supported public input contract is one hierarchical YAML samplesheet.

```yaml
# Generic hierarchical samplesheet template.
# Keep one biological sample entry per sample under `samples:`.
# Omit the `rna:` or `dna:` block if that modality is not present.

library_name: Isa

runtime:
  env_prefix: /path/to/env/prefix

references:
  species: human
  root: /path/to/TrESFlow_References
  ligation_barcode_whitelist: /path/to/TrESFlow_References/ligation_barcode_whitelist.txt
  rna_ref_dir: /path/to/TrESFlow_References/rna/human/star
  dna_ref_dir: /path/to/TrESFlow_References/dna/human/bwa
  dna_blacklist_bed: /path/to/TrESFlow_References/dna/human/hg38-blacklist.v2.bed
  dna_chrom_sizes: /path/to/TrESFlow_References/dna/human/hg38.chrom.sizes
  dna_effective_genome_size: 2913022398

samples:
  sample_id:
    groups:
      group_a:
        rna_sb_barcodes: [AAAA, CCCC]
        dna_sb_barcodes: [AAA, CCC]
      group_b:
        rna_sb_barcodes: [GGGG, TTTT]
        dna_sb_barcodes: [GGG, TTT]

    rna:
      reads:
        i1: /path/to/sample_I1.fq.gz
        r1: /path/to/sample_R1.fq.gz
        r2: /path/to/sample_R2.fq.gz

    dna:
      # Use `single` for legacy 4 nt DNA sample barcodes, or `dual` for explicit
      # 3 nt DNA sample barcodes that differ from RNA barcodes. `reads.i2` is
      # required for single tagmentation and optional for dual tagmentation.
      # Ligation starts are 15,53,91 from i1 in single mode and 41,79,117 from
      # i1 in dual mode.
      tagmentation: dual
      reads:
        i1: /path/to/sample_DNA_I1.fq.gz
        r1: /path/to/sample_DNA_R1.fq.gz
        r2: /path/to/sample_DNA_R2.fq.gz
      mark_barcodes:
        H3K27me3: AGGCTATA
        H3K27ac: GCCTCTAT
```

Notes:

- Omit the `rna:` or `dna:` block when a sample has only one modality. At least one modality block must be present for each sample.
- `runtime.env_prefix`, `references.species`, `references.root`, and `references.ligation_barcode_whitelist` are required. `runtime.tmpdir` is optional and defaults to `--outdir`. Runtime and reference paths are no longer accepted as normal CLI parameters.
- `groups.<group>.sb_barcodes` remains supported for single-tagmentation samples. Use `rna_sb_barcodes` and `dna_sb_barcodes` when RNA and DNA sample barcodes differ; `dna.tagmentation: dual` requires explicit 3 nt `dna_sb_barcodes`.
- `dna.reads.i2` is required for single tagmentation and optional for dual tagmentation.
- DNA ligation tagging uses `reads.i1` starts `15,53,91` for single tagmentation and `reads.i1` starts `41,79,117` for dual tagmentation.
- `references.rna_ref_dir` is required when RNA samples are present and must point directly to the STAR index directory.
- `references.dna_ref_dir`, `references.dna_blacklist_bed`, and `references.dna_effective_genome_size` are required when DNA samples are present.
- `dna.mark_barcodes` is the source of truth for DNA modality barcodes.

Committed examples:

- smoke test: [`assets/samplesheet.example.yaml`](assets/samplesheet.example.yaml)
- canonical real example: [`assets/samplesheet.real.example.yaml`](assets/samplesheet.real.example.yaml)
- generic template: [`assets/samplesheet.template.yaml`](assets/samplesheet.template.yaml)

## Workflow Summary

- [`docs/architecture/implemented_pipeline.md`](docs/architecture/implemented_pipeline.md)

## Outputs

RNA publishes:

- `rna_split_fastqs/`
- `rna_align/`
- `TrES_Stats/`
- `pipeline_info/`

DNA publishes:

- `dna_split_fastqs/`
- `dna_align/`
- `TrES_Stats/`
- `pipeline_info/`

`TrES_Stats/` includes RNA and DNA sequencing-efficiency UpSet PDF plots. Sankey plots, HTML reports, count tables, combined RNA+DNA reports, and sequencing-efficiency warning TSVs are not produced. Optional unavailable BAM-derived categories are skipped with warnings in the process log.

## Runtime Contract

The runtime contract comes from the samplesheet `runtime:` block. `runtime.tmpdir`
is optional; when omitted, `--outdir` is exported as `TMPDIR` for pipeline tasks.
This directory can become very large on real runs. The pipeline creates the directory
if it is missing and fails if it is not writable.

Reference paths are explicit in the samplesheet:

```text
references:
  species: human
  root: /path/to/TrESFlow_References
  ligation_barcode_whitelist: /path/to/TrESFlow_References/ligation_barcode_whitelist.txt
  rna_ref_dir: /path/to/TrESFlow_References/rna/human/star
  dna_ref_dir: /path/to/TrESFlow_References/dna/human/bwa
  dna_blacklist_bed: /path/to/TrESFlow_References/dna/human/hg38-blacklist.v2.bed
  dna_chrom_sizes: /path/to/TrESFlow_References/dna/human/hg38.chrom.sizes
  dna_effective_genome_size: 2913022398
```

`references.rna_ref_dir` is passed directly to STAR as `--genomeDir`. The directory must contain `Genome`, `SA`, `SAindex`, `chrName.txt`, `chrLength.txt`, `chrStart.txt`, `chrNameLength.txt`, and `genomeParameters.txt`.

`references.dna_ref_dir` must contain exactly one complete bwa-mem2 sidecar set. The pipeline infers the prefix from files such as `hg38.fa.0123`, `hg38.fa.amb`, `hg38.fa.ann`, `hg38.fa.bwt.2bit.64`, and `hg38.fa.pac`.

The main remaining runtime CLI parameter is `--max_cpus`, with optional process-bucket CPU overrides for local execution.

Default local CPU budget:

- `--max_cpus 64`
- `--cleanup_work true`, which removes successful task work directories after a successful run.
- `RNA_STARSOLO_ALIGN` defaults to `--rna_starsolo_cpus 16`.
- `ALIGN_DNA` defaults to `--dna_align_cpus 16` and passes that value to bwa-mem2 and samtools.
- `RNA_COVERAGE` and `BAM_COVERAGE_DNA` default to `--coverage_cpus 8`.
- `RNA_FILTERED_BAM`, trim, split, and duplicate-filter helper steps default to `--helper_cpus 4`.
- tagging processes default to `--tagging_cpus 4` and `--tagging_memory '32 GB'`.
- `FQ_TO_SAM` and `MARK_DUPLICATES_DNA` stay at `1` core.
- These are scheduler reservations capped by `--max_cpus`; users can override them with CLI params or Nextflow config.

Work-directory cleanup is intentionally aggressive: `--cleanup_work true` uses Nextflow's successful-run cleanup to remove task work directories after outputs have been published and downstream tasks have completed. This substantially reduces retained `work/` storage, but cleaned tasks are not expected to be usable with `--resume`. Set `--cleanup_work false` when you need the previous resume-friendly behavior for debugging or iterative development.

DNA alignment no longer removes low-count cell barcodes during `ALIGN_DNA`. The aligned BAM still keeps proper-pair mapped, non-blacklisted reads; duplicate removal is represented later by `*_NoDup.bam`, and duplicate status appears in DNA sequencing-efficiency plots as `Unique +`.

Every run writes:

- `${outdir}/pipeline_info/execution_report.html`
- `${outdir}/pipeline_info/execution_timeline.html`
- `${outdir}/pipeline_info/execution_trace.tsv`
- `${outdir}/pipeline_info/flowchart.html`
- `${outdir}/pipeline_info/runtime_contract.tsv`

The active runtime scripts live under [`scripts/core_runtime/`](scripts/core_runtime/). `upstream/source_scripts/` is kept only as provenance for the vendored core code.

## Quick Start

```bash
NXF_OFFLINE=true nextflow run . \
  --samplesheet /path/to/samplesheet.yaml \
  --outdir /path/to/TrESFlow_results \
  --max_cpus 32
```
