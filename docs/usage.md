# TrESFlow: Usage

## Overview

`TrESFlow` supports one public input contract: a single hierarchical YAML samplesheet passed with `--samplesheet`.
There is no CSV input mode in this repository.

The pipeline runs two independent modality branches from that YAML, then builds nf-core QC sidecar outputs, MultiQC, and a self-contained TrESFlow-specific HTML report from explicit metric channels:

- `rna`: sample-barcode tagging, UMI tagging, cell-barcode tagging, trimming, split by SB groups, `FqToSAM`, STARsolo, filtered BAM, bigWigs
- `dna`: sample-barcode tagging, modality tagging, cell-barcode tagging, trimming, split by SB groups and DNA marks, alignment, duplicate marking, NoDup BAM, bigWig

End-of-run reporting is collected under `TrES_Stats`:

- `TrES_Stats/tres_report.html`: offline run overview, per-independent-run retention, sample/mark composition, DNA duplicate complexity, and STARsolo RNA sequencing-saturation report with embedded inline SVG
- `TrES_Stats/read_retention.tsv`, `qc_metrics.tsv`, `barcode_composition.tsv`, and `library_complexity.tsv`: consolidated auditable data
- `TrES_Stats/qc/multiqc/multiqc_report.html`: nf-core MultiQC aggregation of supported logs and QC files
- `TrES_Stats/qc/fastqc/*_fastqc.{html,zip}`: nf-core FastQC reports for raw FASTQs
- `TrES_Stats/qc/samtools/*.flagstat`, `*.stats`, `*.idxstats`, and `*.quickcheck.tsv`: combined Samtools sidecar QC for real BAM outputs

The pipeline report title is the basename of the resolved `--outdir`; its
overview cards retain the samplesheet group names. The standalone assessor uses
the assessed output-directory basename by default and accepts `--title` (or the
backward-compatible `--library-name` alias) as an explicit override.

The samtools sidecars are disabled in `-profile test` because the smoke profile uses mock BAM text files rather than valid BAMs.

DNA alignment no longer filters out low-count cell barcodes during `ALIGN_DNA`. Duplicate-aware DNA outputs are represented by the published `*_MarkedDup.bam` and `*_NoDup.bam` files.

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

TrESFlow supports both Nextflow parser v1 and parser v2 from one source tree
(Nextflow 24.10 or later; parser v2 is the default from Nextflow 26.04).
Paths supplied on the command line or in launch-time config are based on the
launch directory. Therefore relative `--samplesheet`, `--outdir`, and explicit
`--core_scripts_dir` values behave the same whether the pipeline is invoked as
`.` or by an absolute path from another directory. The default output directory
is `<launch-directory>/results`. Repository-owned wrappers, assets, modules,
and the default `scripts/core_runtime` remain based on the pipeline project
directory. Relative paths written inside a samplesheet continue to resolve
from that samplesheet's directory.

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
        mark_barcodes:
          H3K27me3: AGGCTATA
          H3K27ac: GCCTCTAT
      Co2:
        rna_sb_barcodes: [GTCA, TGCA]
        dna_sb_barcodes: [GTC, TGC]
        mark_barcodes:
          H3K27me3: AGGCTATA
          H3K27ac: GCCTCTAT

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
```

### Top-level fields

- `library_name`: run-level library label propagated into RG headers and derived contract files
- `runtime`: required runtime environment and explicit task temporary directory
- `references`: required species label, shared files, and direct RNA/DNA reference paths
- `samples`: biological sample blocks keyed by user-defined sample ID

### `runtime`

- `env_prefix`: environment prefix containing `python3`, `codon`, `cutadapt` (5.2 in the supported environment), `trim_galore`, `STAR`, `samtools`, `bedGraphToBigWig`, `bwa-mem2`, `bamCoverage`, `FastQC`, and `gatk`
- `tmpdir`: optional explicit task temporary directory. If omitted, the pipeline uses `--outdir`. The pipeline creates it if missing and fails if it is not writable.

### `references`

`references.rna_ref_dir` points directly to the STAR index directory. The pipeline passes this exact path to STAR and does not append species, `rna`, or `star`.

`references.dna_ref_dir` points to the directory containing exactly one complete bwa-mem2 sidecar set. The inferred prefix is used for `bwa-mem2 mem`.

`references.dna_effective_genome_size` is required for DNA runs because `DEEPTOOLS_BAMCOVERAGE` passes it to `bamCoverage --effectiveGenomeSize`.

For human references, TrESFlow derives canonical chromosome allowlists once
from STAR's `chrNameLength.txt` and the bwa-mem2 `.ann` dictionary. It accepts a
coherent UCSC convention (`chr1`-`chr22`, `chrX`, `chrY`, `chrM`) or Ensembl
convention (`1`-`22`, `X`, `Y`, and exactly `MT` or `M`). It never renames
contigs or combines conventions. Missing mitochondrial/X/Y anchors, mixed
conventions, or disagreement between the DNA `.ann` dictionary and an explicit
`dna_chrom_sizes` file are fatal configuration errors.

The resulting `*_canonical_chromosomes.txt` and
`*_canonical_chromosomes.chrom.sizes` files are recorded in
`pipeline_info/derived_contract/`. RNA's canonical filtered BAM feeds QC,
normal-compression publication, and stranded/unstranded coverage. DNA keeps
duplicate marking unchanged, filters the normalized MarkedDup BAM afterward,
and derives the final canonical NoDup BAM and BigWig from it.

### `samples.<sample_id>.groups`

`groups` is the source of truth for biological sample-barcode grouping.

- each group key is the biological label that will appear in split outputs
- `rna_sb_barcodes` and `dna_sb_barcodes` are modality-specific sample barcodes assigned to that logical group.
- A group participates in RNA when it has `rna_sb_barcodes`; it participates in DNA when it has both `dna_sb_barcodes` and `mark_barcodes`. RNA-only and DNA-only groups may coexist under one sample.
- Legacy `sb_barcodes` remains supported for single-tagmentation samples. In `dna.tagmentation: dual`, DNA requires explicit 3 nt `dna_sb_barcodes`; they are not derived from RNA barcodes.
- `mark_barcodes` belongs inside each DNA-participating group and maps mark labels to modality barcodes for that group.
- sample barcodes must be unique within a sample block

### `samples.<sample_id>.rna`

The RNA block is optional, but if present it must contain:

- `reads.i1`
- `reads.r1`
- `reads.r2`

### `samples.<sample_id>.dna`

The DNA block is optional, but if present it must contain:

- `tagmentation`: either `single` or `dual` (required; there is no implicit default)
- `reads.i1`
- `reads.r1`
- `reads.r2`

`reads.i2` is required for `dna.tagmentation: single` and optional for `dna.tagmentation: dual`.

DNA ligation tagging uses the same `Tag_Lig3` correction and output format for both modes, but with mode-specific barcode-source reads and start positions:

- `single`: ligation source `reads.i1`, L1/L2/L3 starts `15,53,91`
- `dual`: ligation source `reads.i1`, L1/L2/L3 starts `41,79,117`

`groups.<group>.mark_barcodes` maps biological mark labels to DNA modality
barcodes for that group. A group participates in RNA when it defines
`rna_sb_barcodes` (with legacy `sb_barcodes` retained as a single-tag fallback),
and participates in DNA when it defines both `dna_sb_barcodes` and
`mark_barcodes`. RNA-only and DNA-only groups may coexist under one sample;
the sample-level `rna.reads` and `dna.reads` are each processed once.

## Derived Internal Contract

The parser in [`lib/SamplesheetParser.groovy`](../lib/SamplesheetParser.groovy) turns the hierarchical YAML into modality-specific work rows and writes helper contract files under:

```text
<outdir>/pipeline_info/derived_contract/
```

These derived files include:

- `sb_group_map.tsv`
- `dna_mo_map.tsv` when DNA is present
- per-sample DNA modality whitelist files

`dna_mo_map.tsv` retains the four columns `sample`, `sb_group`, `mark`, and
`mo_bc`; mappings are group-specific, so one MO barcode may identify different
marks in different groups. The per-sample DNA modality whitelist is the union
of MO barcodes across that sample's DNA-participating groups.

This keeps the public input contract user-friendly while preserving the split and alignment interfaces used by the current core modules.

## Parameters

The main public parameters are:

- `--samplesheet`
- `--outdir`
- `--publish_split_fastqs`
- `--max_cpus`
- `--cleanup_work`
- `--rna_starsolo_cpus`
- `--dna_align_cpus`
- `--coverage_cpus`
- `--helper_cpus`
- `--tagging_cpus`
- `--tagging_memory`
- `--filter_dual_tag_artifacts`
- `--aviti_optical_duplicate_distance`

Deprecated runtime/reference CLI parameters now fail with a hard error.

`--filter_dual_tag_artifacts` defaults to `true`. For dual-tagmentation DNA
only, TrESFlow searches both mates for 48 audited, oriented linker signatures
after cell-barcode tagging and Trim Galore. Trimming occurs first so a valid
genomic prefix can be salvaged before downstream linker read-through is
examined. Each signature is an exact 23-mer; a full match remaining anywhere
in either trimmed mate discards the pair, while retained FASTQ records pass
through unchanged relative to the Trim Galore outputs. Single-tagmentation DNA
and all RNA samples bypass this process, as do dual samples when the parameter
is `false`.
The version-controlled signature asset has SHA-256
`67c6f1789ef5e36492562203ac38fc13fa901058047ed2bd37b304d85a30ae0f`;
the supplied provenance audit records zero exact matches in hg38, mm10, and
mm39. The filter has no QC pass/fail threshold. Trim Galore remains responsible
for ordinary adapter/quality trimming, quality threshold 10, and the 20-base
minimum before exact residual-signature filtering begins.

The artifact summary's `input_pairs` value is therefore the number of paired
reads surviving Trim Galore, not the raw or pre-trimming pair count.

The published `*.barcode_gates.tsv` files provide exact cumulative same-pair
barcode populations among reads that reach splitting. They are computed from
the routing decisions already recorded during tagging, not by rematching
barcodes. The companion `*.barcode_composition.tsv` files give one exact
category per configured sample-barcode sequence with its group/label, or per
configured DNA mark, plus `NoMatch` even when zero, with explicit denominators.
Marginal barcode percentages remain raw-input diagnostics and
must not be interpreted as sequential retention stages. In an enabled dual-tag
DNA run, `split_input_pairs` and every barcode-gate denominator begin after the
artifact filter; otherwise they begin after paired trimming.

`--aviti_optical_duplicate_distance` is a non-negative integer and defaults to
`10`. It is the empirical AVITI coordinate-distance threshold passed to Picard
MarkDuplicates for spatial/optical duplicate classification. The value was
calibrated on AVITI 500 data, is expressed in AVITI coordinate units, and is
not an Illumina pixel-distance setting or a universal Element-defined value.
Override it only for an independently calibrated dataset.

DNA duplicate grouping remains cell-aware through the corrected 24-base `CB`
tag. AVITI read names preserve `lane:tile:x:y`; TrESFlow parses those names with
an AVITI-specific Picard regex and writes a compact lane-level DNA read-group
namespace (`L1`, `L2`, ...). All lane read groups for a logical library have the
same `LB`, so ordinary genomic/PCR duplicate families can span lanes for the same
cell while Picard's physical comparison cannot collide coordinates from two
different lanes. Cell identity is not encoded in DNA `RG`; it remains in `CB`.
RNA read-group behavior is unchanged.

The standalone assessor uses the same repository-owned parser, validation,
normalized model, plots, and HTML renderer as the pipeline:

```bash
python3 bin/assess_tresflow_run.py /path/to/completed/output \
  --output-dir /path/to/qc-assessment
```

It is read-only with respect to the assessed output. Legacy runs without exact
gate sidecars retain only provable coarse stages; exact gates or composition are
never reconstructed from marginal percentages.

Unlike the pipeline report process, the standalone assessor also writes its
normalized JSON model and standalone SVG copies by default. Use `--no-json` and
`--no-standalone-figures` to request the pipeline's compact five-file contract.

`--publish_split_fastqs` defaults to `false`. Internal plain split FASTQs still
feed RNA and DNA computation, but the publication-only `pigz` tasks and the
top-level `rna_split_fastqs/` and `dna_split_fastqs/` directories are omitted.
Set `--publish_split_fastqs` to publish both modalities' split FASTQs, when
present, in their existing gzip-compressed directories and filenames.

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
