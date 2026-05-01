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
git clone git@github.com:AAnnan/TrESFlow.git
cd TrESFlow
```

Install codon in your env:
```bash
./scripts/install_codon_0.16.3.sh --prefix /path/to/env/prefix (for ex:/home/ahrmad/micromamba/envs/tres)
```

## Samplesheet Contract

The only supported public input contract is one hierarchical YAML samplesheet.

```yaml
library_name: Isa

runtime:
  env_prefix: /home/annan/micromamba/envs/tres
  tmpdir: /mnt/dataFast/ahrmad/tmp/TrESFlow_Isa

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
        sb_barcodes: [CAGT, ACGT, GCTA, CGTA]
      Co2:
        sb_barcodes: [GTCA, TGCA, CTGA, TCGA]

    rna:
      reads:
        i1: test_realdata/day15_I1.fq.gz
        r1: test_realdata/day15_R1.fq.gz
        r2: test_realdata/day15_R2.fq.gz

    dna:
      reads:
        i1: test_realdata/day15_DNA_I1.fq.gz
        i2: test_realdata/day15_DNA_I2.fq.gz
        r1: test_realdata/day15_DNA_R1.fq.gz
        r2: test_realdata/day15_DNA_R2.fq.gz
      mark_barcodes:
        H3K27me3: AGGCTATA
        H3K27ac: GCCTCTAT
```

Notes:

- Omit the `rna:` or `dna:` block when a sample has only one modality. At least one modality block must be present for each sample.
- `runtime.env_prefix`, `runtime.tmpdir`, `references.species`, `references.root`, and `references.ligation_barcode_whitelist` are required. Runtime and reference paths are no longer accepted as normal CLI parameters.
- `groups.<group>.sb_barcodes` remains supported for single-tagmentation samples. Use `rna_sb_barcodes` and `dna_sb_barcodes` when RNA and DNA sample barcodes differ; `dna.tagmentation: dual` requires explicit 3 nt `dna_sb_barcodes`.
- `references.rna_ref_dir` is required when RNA samples are present and must point directly to the STAR index directory.
- `references.dna_ref_dir`, `references.dna_blacklist_bed`, and `references.dna_effective_genome_size` are required when DNA samples are present.
- `dna.mark_barcodes` is the source of truth for DNA modality barcodes.

Committed examples:

- smoke test: [`assets/samplesheet.example.yaml`](assets/samplesheet.example.yaml)
- canonical real example: [`assets/samplesheet.real.example.yaml`](assets/samplesheet.real.example.yaml)
- generic template: [`assets/samplesheet.template.yaml`](assets/samplesheet.template.yaml)

## Workflow Summary

RNA core:

1. `TAG_RNA_SAMPLE_BARCODE`
2. `TAG_RNA_UMI`
3. `TAG_RNA_CELL_BARCODE`
4. `TRIM_RNA_FASTQS`
5. `SPLIT_RNA_READS`
6. `FQ_TO_SAM`
7. `RNA_STARSOLO_ALIGN`
8. `RNA_FILTERED_BAM`
9. `RNA_COVERAGE`

DNA core:

1. `TAG_DNA_SAMPLE_BARCODE`
2. `TAG_DNA_MODALITY_BARCODE`
3. `TAG_DNA_CELL_BARCODE`
4. `TRIM_DNA_FASTQS`
5. `SPLIT_DNA_READS`
6. `ALIGN_DNA`
7. `MARK_DUPLICATES_DNA`
8. `SPLIT_DUPLICATES_DNA`
9. `BAM_COVERAGE_DNA`

Architecture/DAG:

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

## Runtime Contract

The runtime contract comes from the samplesheet `runtime:` block. `runtime.tmpdir`
is exported as `TMPDIR` for pipeline tasks and can become very large on real runs.
The pipeline creates the directory if it is missing and fails if it is not writable.

Reference paths are explicit in the samplesheet:

```yaml
references:
  species: human
  root: /mnt/dataFast/ahrmad/TrESFlow_References
  ligation_barcode_whitelist: /mnt/dataFast/ahrmad/TrESFlow_References/ligation_barcode_whitelist.txt
  rna_ref_dir: /mnt/dataFast/ahrmad/TrESFlow_References/rna/human/star
  dna_ref_dir: /mnt/dataFast/ahrmad/TrESFlow_References/dna/human/bwa
  dna_blacklist_bed: /mnt/dataFast/ahrmad/TrESFlow_References/dna/human/hg38-blacklist.v2.bed
  dna_chrom_sizes: /mnt/dataFast/ahrmad/TrESFlow_References/dna/human/hg38.chrom.sizes
  dna_effective_genome_size: 2913022398
```

`references.rna_ref_dir` is passed directly to STAR as `--genomeDir`. The directory must contain `Genome`, `SA`, `SAindex`, `chrName.txt`, `chrLength.txt`, `chrStart.txt`, `chrNameLength.txt`, and `genomeParameters.txt`.

`references.dna_ref_dir` must contain exactly one complete bwa-mem2 sidecar set. The pipeline infers the prefix from files such as `hg38.fa.0123`, `hg38.fa.amb`, `hg38.fa.ann`, `hg38.fa.bwt.2bit.64`, and `hg38.fa.pac`.

The main remaining runtime CLI parameter is `--max_cpus`.

Default local CPU budget:

- `--max_cpus 64`
- `RNA_STARSOLO_ALIGN` reserves up to `48` cores.
- `RNA_COVERAGE` and `BAM_COVERAGE_DNA` reserve up to `32` cores.
- `ALIGN_DNA` reserves up to `48` cores and passes that value to bwa-mem2 and samtools.
- `RNA_FILTERED_BAM`, trim, split, and duplicate-filter helper steps reserve up to `16` cores.
- tagging processes default to `6` cores and `64 GB` memory.
- `FQ_TO_SAM` and `MARK_DUPLICATES_DNA` stay at `1` core.
- These are scheduler reservations derived from `--max_cpus`; Nextflow still prevents local tasks from exceeding the configured CPU budget.

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
  --samplesheet /mnt/dataFast/ahrmad/TEST_NF/isa_multiome.yaml \
  --outdir /mnt/dataFast/ahrmad/TEST_NF/TrESFlow_Isa \
  --max_cpus 32
```
