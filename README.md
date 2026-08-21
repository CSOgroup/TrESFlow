# TrESFlow

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)

TrESFlow is a Nextflow DSL2 pipeline for preprocessing TrES-seq RNA and DNA sequencing data, including barcode assignment, read splitting, alignment, duplicate handling, coverage tracks, and quality control.

## Installation

Create and activate a conda/mamba/micromamba environment with the required tools:

```bash
micromamba create -n tres
micromamba activate tres

micromamba install \
  pandas polars ipython pysam pybedtools numpy matplotlib seaborn scipy \
  pyarrow upsetplot anndata scanpy matplotlib-venn leidenalg scikit-learn \
  snapatac2 screen samtools bwa-mem2 star fastqc multiqc trim-galore \
  deeptools parallel ucsc-bedGraphToBigWig nextflow git gatk4
```

Clone TrESFlow:

```bash
git clone git@github.com:CSOgroup/TrESFlow.git
cd TrESFlow
```

Install Codon in the same environment:

```bash
./scripts/install_codon_0.16.3.sh --prefix /path/to/env/prefix
```

## Inputs

TrESFlow uses one hierarchical YAML samplesheet describing the runtime environment, references, biological groups, barcodes, and FASTQ inputs.

```yaml
library_name: Isa

runtime:
  env_prefix: /path/to/env/prefix
  # tmpdir: /path/to/large/tmp

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
        mark_barcodes:
          H3K27me3: AGGCTATA
          H3K27ac: GCCTCTAT

      group_b:
        rna_sb_barcodes: [GGGG, TTTT]
        dna_sb_barcodes: [GGG, TTT]
        mark_barcodes:
          H3K27me3: AGGCTATA
          H3K9me3: GCCTCTAT

    rna:
      reads:
        i1: /path/to/sample_RNA_I1.fastq.gz
        r1: /path/to/sample_RNA_R1.fastq.gz
        r2: /path/to/sample_RNA_R2.fastq.gz

    dna:
      # Both `single` and `dual` tagmentation modes are supported.
      #
      # Single tagmentation uses 4 nt DNA sample barcodes.
      # Dual tagmentation uses 3 nt DNA sample barcodes.
      #
      # `reads.i2` is required for single tagmentation, not for dual.
      #
      # DNA ligation barcodes are read from i1 at positions:
      #   single: 15,53,91
      #   dual:   41,79,117
      tagmentation: dual
      reads:
        i1: /path/to/sample_DNA_I1.fastq.gz
        r1: /path/to/sample_DNA_R1.fastq.gz
        r2: /path/to/sample_DNA_R2.fastq.gz
```

### Samplesheet notes

- A sample may contain RNA, DNA, or both.
- FASTQs are defined once at sample level. Groups define how reads are assigned downstream.
- Groups may independently participate in RNA and/or DNA.
- RNA groups use `rna_sb_barcodes`.
- DNA groups use `dna_sb_barcodes` together with `mark_barcodes`.
- The same DNA modality barcode may represent different marks in different groups.
- `dna.reads.i2` is required for `single` tagmentation and optional for `dual`.
- Every `reads.i1`, `reads.i2`, `reads.r1`, and `reads.r2` value may be one path, a comma-separated scalar, or a YAML sequence. Lists are only for ordered technical FASTQ chunks of the same biological library, chemistry, tagmentation mode, and barcode layout.
- Entry 1 in each read role is one synchronized technical read set, entry 2 is the next set, and so on. Required roles must have equal lengths; TrESFlow does not scan directories, expand globs, or infer mates from filenames.
- Use a YAML sequence when a filename itself contains a comma. Sequence order is preserved. Relative entries resolve from the samplesheet directory.
- `runtime.tmpdir` is optional and defaults to `--outdir`.
- `references.rna_ref_dir` must point to a STAR index when RNA is present.
- DNA samples require a bwa-mem2 reference, blacklist, chromosome sizes, and effective genome size.

For example, two RNA run chunks can be written as comma-separated scalars:

```yaml
rna:
  reads:
    i1: /run1/I1.fastq.gz, /run2/I1.fastq.gz
    r1: /run1/R1.fastq.gz, /run2/R1.fastq.gz
    r2: /run1/R2.fastq.gz, /run2/R2.fastq.gz
```

The equivalent YAML-sequence form is the escape hatch for paths containing commas. TrESFlow stages all sources collision-safely, streams them in configured order into one set of tagged outputs, and keeps downstream filenames and reports at one logical-library level. Every explicitly supplied raw FASTQ reaches FastQC once; an omitted dual-DNA `i2` uses the internal `i1` fallback without a duplicate FastQC run.

An example samplesheet is available above and a template in: [`assets/samplesheet.template.yaml`](assets/samplesheet.template.yaml)

## Running TrESFlow

A typical run is:

```bash
NXF_OFFLINE=true nextflow run . \
  --samplesheet /path/to/samplesheet.yaml \
  --outdir /path/to/results \
  --max_cpus 80
```

Useful options include:

- `--max_cpus`: maximum CPU budget available to the pipeline.
- `--cleanup_work false`: retain successful Nextflow work directories when debugging or when `-resume` is important.
- `--publish_split_fastqs true`: also publish gzip-compressed per-group split FASTQs.
- `--filter_dual_tag_artifacts false`: disable the dual-tagmentation DNA artifact filter.

`NXF_OFFLINE=true` prevents loading optional remote nf-core configuration and is suitable for normal local TrESFlow runs.

## Outputs

The main output structure is:

```text
<outdir>/
├── rna_align/
│   ├── STARsolo outputs
│   ├── *.filtered_cells.bam
│   └── *.bw
│
├── dna_align/
│   ├── aligned BAMs
│   ├── *_MarkedDup.bam
│   ├── *_NoDup.bam
│   └── *_NoDup.bw
│
├── TrES_Stats/
│   ├── tres_report.html
│   ├── read_retention.tsv
│   ├── qc_metrics.tsv
│   ├── barcode_composition.tsv
│   ├── library_complexity.tsv
│   └── qc/
│       ├── fastqc/
│       ├── samtools/
│       └── multiqc/
│
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.tsv
    ├── flowchart.html
    └── derived_contract/
```

The main QC report is:

```text
TrES_Stats/tres_report.html
```

It summarizes read retention, barcode composition, DNA duplicate/library-complexity metrics, and RNA sequencing metrics.

MultiQC is written to:

```text
TrES_Stats/qc/multiqc/multiqc_report.html
```

When `--publish_split_fastqs true` is used, the pipeline additionally creates:

```text
rna_split_fastqs/
dna_split_fastqs/
```

## Important information

### Work directory cleanup

`--cleanup_work true` is the default. After a successful run, completed task work directories are removed to save disk space.

Use:

```bash
--cleanup_work false
```

when you need to preserve the work directory for debugging or reliable reuse with `-resume`.

### Temporary disk space

`runtime.tmpdir` defaults to `--outdir`. Large datasets can require substantial temporary disk space, so set `runtime.tmpdir` explicitly when a larger or faster filesystem is available.

### Canonical chromosomes

TrESFlow derives canonical chromosomes directly from the supplied references. UCSC and Ensembl chromosome naming conventions are supported.

Final canonical RNA and DNA outputs exclude alternative loci, patches, random/unplaced contigs, decoys, and other noncanonical sequences.

### DNA duplicates

DNA duplicate marking uses the corrected cell barcode and AVITI read coordinates. DNA `RG` and `PU` identify the full physical unit (`instrument:run:flowcell:L<lane>`), while all units retain the same logical `LB`. This allows PCR/library duplicate families to span units but keeps optical comparisons unit-local. Unsupported identifier characters fail instead of being normalized into a colliding ID. The default optical-duplicate distance is:

```text
--aviti_optical_duplicate_distance 10
```

Final `*_NoDup.bam` files and DNA BigWig tracks exclude reads marked as duplicates.

### Dual-tagmentation artifact filtering

Dual-tagmentation DNA libraries are filtered after trimming for known residual linker signatures. If either mate contains one of these signatures, the pair is discarded.

This filter:

- is enabled by default;
- applies only to `dual` DNA libraries;
- does not affect RNA or `single` DNA libraries.

## Further documentation

More detailed documentation is available in:

- [Usage documentation](docs/usage.md)
- [Pipeline architecture](docs/architecture/implemented_pipeline.md)
