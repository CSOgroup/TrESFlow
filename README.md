# TrESFlow

TrESFlow is a Nextflow DSL2 pipeline for the implemented TrES core workflow in this repo: RNA through the repo-owned alignment/filtered-BAM/coverage path and DNA through `BAM_COVERAGE_DNA`.

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

## Samplesheet Contract

The only supported public input contract is one hierarchical YAML samplesheet.

```yaml
library_name: REALDATATESTLIB

resources:
  ligation_barcode_whitelist: test_realdata/ligation_barcode_whitelist.txt
  rna_ref_base_dir: /path/to/reference_base_dir
  rna_align_species: human
  dna_bwa_reference: /path/to/bwa_index_prefix
  dna_blacklist_bed: /path/to/blacklist.bed
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

- `resources` holds shared run resources for the whole run.
- Each biological sample appears once under `samples:`.
- Omit the `rna:` or `dna:` block when a sample has only one modality.
- `groups.<group>.sb_barcodes` is the source of truth for sample-barcode grouping.
- `dna.mark_barcodes` is the source of truth for DNA modality barcodes.
- `sb_group_map.tsv`, `dna_mo_map.tsv`, and per-sample DNA modality whitelist files are derived internally under `${outdir}/pipeline_info/derived_contract/`.
- At least one modality block must be present for each sample.

Committed examples:

- smoke test: [`assets/samplesheet.example.yaml`](assets/samplesheet.example.yaml)
- canonical real example: [`assets/samplesheet.real.example.yaml`](assets/samplesheet.real.example.yaml)
- generic template: [`assets/samplesheet.template.yaml`](assets/samplesheet.template.yaml)

## Runtime Contract

The primary runtime contract is one environment root:

- `--runtime_env_prefix`

On this server it defaults to:

- `/home/annan/micromamba/envs/tres`

The pipeline derives these standard executables from `${runtime_env_prefix}/bin`:

- `python3`
- `trim_galore`
- `STAR`
- `samtools`
- `bedGraphToBigWig`
- `bwa-mem2`
- `bamCoverage`
- `gatk`
- `codon`

Seq `0.11.3` is expected under `${runtime_env_prefix}/lib/codon/plugins/seq`.

`nextflow` itself should be launched from the intended environment. The workflow does not choose its own launcher binary.

Shared run resources come from the samplesheet `resources:` block. CLI params remain available only as explicit overrides:

- `--ligation_barcode_whitelist`
- `--rna_ref_base_dir`
- `--rna_align_species`
- `--dna_bwa_reference`
- `--dna_blacklist_bed`
- `--dna_effective_genome_size`
- `--max_cpus`

Default local CPU budget:

- `--max_cpus 64`
- `RNA_STARSOLO_ALIGN` and `RNA_COVERAGE` reserve up to `24` cores each.
- `RNA_FILTERED_BAM`, trim, split, duplicate-filter, and DNA coverage helper steps reserve up to `8` cores each.
- `ALIGN_DNA` reserves up to `32` cores.
- tagging, `FQ_TO_SAM`, and `MARK_DUPLICATES_DNA` stay at `1` core.
- These are scheduler reservations. `ALIGN_DNA` still wraps an upstream script with its own internal thread behavior, so its Nextflow CPU value primarily controls local concurrency.

Every run writes:

- `${outdir}/pipeline_info/execution_report.html`
- `${outdir}/pipeline_info/execution_timeline.html`
- `${outdir}/pipeline_info/execution_trace.tsv`
- `${outdir}/pipeline_info/flowchart.html`
- `${outdir}/pipeline_info/runtime_contract.tsv`

The active runtime scripts live under [`scripts/core_runtime/`](scripts/core_runtime/). `upstream/source_scripts/` is kept only as provenance for the vendored core code.

## Quick Start

Smoke test:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Canonical real example:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.real.example.yaml \
  --outdir results/test_real
```

## Outputs

RNA publishes:

- `tagging/`
- `split/`
- `usam/`
- `align/`
- `pipeline_info/`

DNA publishes:

- `dna_tagging/`
- `dna_split/`
- `dna_align/`
- `dna_dedup/`
- `dna_nodup/`
- `dna_coverage/`
- `pipeline_info/`

## Troubleshooting

- Missing RNA alignment resources: set `resources.rna_ref_base_dir` and `resources.rna_align_species` in the YAML, or use the matching CLI overrides.
- Missing DNA alignment resources: set `resources.dna_bwa_reference`, `resources.dna_blacklist_bed`, and `resources.dna_effective_genome_size` in the YAML, or use the matching CLI overrides.
- Codon/Seq preflight failures: confirm `${runtime_env_prefix}/bin/codon` is `0.16.3` and `${runtime_env_prefix}/lib/codon/plugins/seq` is `0.11.3`. Use `scripts/install_codon_0.16.3.sh` to install both into the environment prefix.
- Long `bamCoverage` tasks on this server are treated as a runtime/performance characteristic, not a pipeline logic blocker, unless they fail with a distinct configuration or execution error.
