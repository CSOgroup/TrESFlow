# TrESFlow

TrESFlow is a Nextflow DSL2 pipeline for the implemented TrES core workflow in this repo: RNA through `ALIGN_RNA` and DNA through `BAM_COVERAGE_DNA`.

## Workflow Summary

RNA core:

1. `TAG_RNA_SAMPLE_BARCODE`
2. `TAG_RNA_UMI`
3. `TAG_RNA_CELL_BARCODE`
4. `TRIM_RNA_FASTQS`
5. `SPLIT_RNA_READS`
6. `FQ_TO_SAM`
7. `ALIGN_RNA`

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
- `groups.<group>.sb_barcodes` is the source of truth for sample-barcode grouping.
- `dna.mark_barcodes` is the source of truth for DNA modality barcodes.
- `sb_group_map.tsv`, `dna_mo_map.tsv`, and per-sample DNA modality whitelist files are derived internally under `${outdir}/pipeline_info/derived_contract/`.
- At least one modality block must be present for each sample.

Committed examples:

- smoke test: [`assets/samplesheet.example.yaml`](assets/samplesheet.example.yaml)
- real RNA template: [`assets/samplesheet.real_rna.template.yaml`](assets/samplesheet.real_rna.template.yaml)
- real RNA example: [`assets/test_realdata/samplesheet.real_rna.yaml`](assets/test_realdata/samplesheet.real_rna.yaml)
- real DNA example: [`assets/samplesheet.dna.RealDATAexample.yaml`](assets/samplesheet.dna.RealDATAexample.yaml)

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

Explicit exceptions that remain outside that env:

- `--runtime_codon`
- `--codon_home`

On this server those default to:

- Codon `0.16.3`: `/home/annan/.codon/bin/codon`
- Seq `0.11.3`: `/home/annan/.codon/lib/codon/plugins/seq`

`nextflow` itself should be launched from the intended environment. The workflow does not choose its own launcher binary.

Shared run resources come from the samplesheet `resources:` block. CLI params remain available only as explicit overrides:

- `--ligation_barcode_whitelist`
- `--rna_ref_base_dir`
- `--rna_align_species`
- `--dna_bwa_reference`
- `--dna_blacklist_bed`
- `--dna_effective_genome_size`
- `--max_cpus`

Every run writes:

- `${outdir}/pipeline_info/execution_report.html`
- `${outdir}/pipeline_info/execution_timeline.html`
- `${outdir}/pipeline_info/execution_trace.tsv`
- `${outdir}/pipeline_info/flowchart.html`
- `${outdir}/pipeline_info/runtime_contract.tsv`

## Quick Start

Smoke test:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Real RNA:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet assets/test_realdata/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna
```

Real DNA:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.dna.RealDATAexample.yaml \
  --outdir results/test_dna_real
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
- Codon/Seq preflight failures: confirm `/home/annan/.codon/bin/codon` is `0.16.3` and the Seq plugin under `/home/annan/.codon/lib/codon/plugins/seq` is `0.11.3`, or override `--runtime_codon` and `--codon_home`.
- Long `bamCoverage` tasks on this server are treated as a runtime/performance characteristic, not a pipeline logic blocker, unless they fail with a distinct configuration or execution error.
