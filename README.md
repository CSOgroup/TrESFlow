# TrESFlow

TrESFlow is a Nextflow DSL2 pipeline that owns the implemented core runtime in
this repo.

- Core workflow:
  - RNA through `ALIGN_RNA`
  - DNA through `BAM_COVERAGE_DNA`
- Optional downstream:
  - shared staging
  - one optional `sc_process.py` call

Architecture/DAG:
- [`docs/architecture/implemented_pipeline.md`](docs/architecture/implemented_pipeline.md)

## Workflow Boundaries

### RNA core

1. `TAG_RNA_SAMPLE_BARCODE`
2. `TAG_RNA_UMI`
3. `TAG_RNA_CELL_BARCODE`
4. `TRIM_RNA_FASTQS`
5. `SPLIT_RNA_READS`
6. `FQ_TO_SAM`
7. `ALIGN_RNA`

### DNA core

1. `TAG_DNA_SAMPLE_BARCODE`
2. `TAG_DNA_MODALITY_BARCODE`
3. `TAG_DNA_CELL_BARCODE`
4. `TRIM_DNA_FASTQS`
5. `SPLIT_DNA_READS`
6. `ALIGN_DNA`
7. `MARK_DUPLICATES_DNA`
8. `SPLIT_DUPLICATES_DNA`
9. `BAM_COVERAGE_DNA`

### Optional downstream

- `STAGE_SC_PROCESS_INPUTS`
- `RUN_SC_PROCESS`

The core workflow is complete without `RUN_SC_PROCESS`. The workflow design
allows only one `sc_process.py` call.

## Validation Policy

Real-data runs are the primary validation target for this repo.

The only retained non-real path is the minimal `-profile test` smoke test,
because [`AGENTS.md`](AGENTS.md) requires the repo to keep an end-to-end
`-profile test` path with a minimal dataset.

Current validation focus:

- real RNA
- real DNA
- optional real multiome staging and downstream analysis

## Quick Start

Minimal required smoke test:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Real RNA validation:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet assets/test_realdata/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna
```

Real DNA validation:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.dna.RealDATAexample.yaml \
  --outdir results/test_dna_real
```

Optional real shared staging:

```bash
nextflow run . \
  --samplesheet assets/test_realdata/samplesheet.real_multiome.yaml \
  --outdir results/test_real_multiome_stage \
  --stage_sc_process_inputs true
```

Optional real downstream analysis:

```bash
nextflow run . \
  --samplesheet assets/test_realdata/samplesheet.real_multiome.yaml \
  --outdir results/test_real_multiome_sc_process \
  --run_sc_process true
```

## Runtime Contract

Default executable bindings on this server:

- `python3`: `/home/annan/micromamba/envs/tres/bin/python3`
- `trim_galore`: `/home/annan/micromamba/envs/tres/bin/trim_galore`
- `STAR`: `/home/annan/micromamba/envs/tres/bin/STAR`
- `samtools`: `/home/annan/micromamba/envs/tres/bin/samtools`
- `bedGraphToBigWig`: `/home/annan/micromamba/envs/tres/bin/bedGraphToBigWig`
- `bwa-mem2`: `/home/annan/micromamba/envs/tres/bin/bwa-mem2`
- `bamCoverage`: `/home/annan/micromamba/envs/tres/bin/bamCoverage`
- `gatk`: `/home/annan/micromamba/envs/tres/bin/gatk`

Current explicit exception:

- Codon `0.16.3`: `/home/annan/.codon/bin/codon`
- Seq `0.11.3`: `/home/annan/.codon/lib/codon/plugins/seq`

Every run writes:

- `${outdir}/pipeline_info/execution_report.html`
- `${outdir}/pipeline_info/execution_timeline.html`
- `${outdir}/pipeline_info/execution_trace.tsv`
- `${outdir}/pipeline_info/flowchart.html`
- `${outdir}/pipeline_info/runtime_contract.tsv`

The implemented core workflow runs from [`scripts/core_runtime/`](scripts/core_runtime/).
[`upstream/source_scripts/`](upstream/source_scripts/) remains for provenance
and for the optional downstream `sc_process.py` entrypoint only.

## Inputs and Params

### Required for every run

- `--samplesheet`
- `--outdir`

Pipeline input is a single YAML samplesheet under `params.samplesheet`.

### Samplesheet contract

The supported public contract is a hierarchical YAML where each biological
sample appears once, with nested group definitions and nested RNA/DNA modality
blocks.

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

User-supplied top-level keys:

- `library_name`
- `resources`
- `samples`

`resources` holds shared run inputs such as:

- `ligation_barcode_whitelist`
- `rna_ref_base_dir`
- `rna_align_species`
- `dna_bwa_reference`
- `dna_blacklist_bed`
- `dna_effective_genome_size`

Per sample:

- `groups.<group>.sb_barcodes`
- optional `rna.reads`
- optional `dna.reads`
- optional `dna.mark_barcodes`

At least one modality block must be present for each sample.

### Internally derived contract files

The pipeline now derives the implementation-facing files it still needs under:

- `${outdir}/pipeline_info/derived_contract/sb_group_map.tsv`
- `${outdir}/pipeline_info/derived_contract/dna_mo_map.tsv`
- `${outdir}/pipeline_info/derived_contract/dna_modality_whitelists/<sample>.txt`

These are derived from the hierarchical YAML and are not part of the public
input contract.

Users no longer supply separate:

- sample-barcode group TSV files
- DNA modality-map TSV files
- DNA modality-barcode whitelist files

The shared sample-barcode groups in the YAML are now the source of truth for:

- grouped RNA splitting
- grouped DNA splitting
- the effective DNA sample-barcode whitelist

The DNA `mark_barcodes` mapping in the YAML is now the source of truth for:

- the internally derived `dna_mo_map.tsv`
- per-sample DNA modality whitelists

### Required for real RNA runs

Required in the run contract:

- `resources.rna_ref_base_dir` in the samplesheet, or `--rna_ref_base_dir` as a fallback
- `resources.rna_align_species` in the samplesheet, or `--rna_align_species` as a fallback

Reference layout:

- human: `<ref_base_dir>/GRCh38_TrES/star` and `<ref_base_dir>/hg38.chrom.sizes`
- mouse: `<ref_base_dir>/GRCm39_TrES/star` and `<ref_base_dir>/mm39.chrom.sizes`

### Required for real DNA runs

Required in the run contract:

- `resources.dna_bwa_reference` in the samplesheet, or `--dna_bwa_reference` as a fallback
- `resources.dna_blacklist_bed` in the samplesheet, or `--dna_blacklist_bed` as a fallback
- `resources.dna_effective_genome_size` in the samplesheet, or `--dna_effective_genome_size` as a fallback

`--dna_bwa_reference` is a bwa-mem2 index prefix. The following sidecars must
exist: `.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`.

### Defaulted barcode resources and tag settings

The ligation barcode whitelist is now part of the preferred samplesheet
contract under `resources.ligation_barcode_whitelist`.

If it is omitted from the YAML, the fallback is:

- `--ligation_barcode_whitelist`

The default server config points it at:

- `assets/test_realdata/ligation_barcode_whitelist.txt`

Codon barcode-tagging details are also kept out of the samplesheet by default.
They are defined in `params.barcode_defaults` in [`nextflow.config`](nextflow.config)
and cover:

- RNA sample barcode settings
- RNA UMI settings
- RNA cell barcode settings
- DNA sample barcode settings
- DNA modality barcode settings
- DNA cell barcode settings

### Optional downstream controls

- `--stage_sc_process_inputs`
- `--run_sc_process`
- `--runtime_snap_data_dir`
- `--optional_sc_process_script`

`RUN_SC_PROCESS` binds `SNAP_DATA_DIR` explicitly. By default on this server it
uses `/home/annan/.cache/snapatac2`, so `snap.genome.hg38` and
`snap.genome.mm39` resolve from the local SnapATAC cache.

### Runtime and layout overrides

- `--core_scripts_dir`
- `--runtime_bin_dir`
- individual `--runtime_*` params
- `--runtime_codon`
- `--codon_home`
- `--max_cpus`

`--max_cpus` defaults to `40`.

The samplesheet is the primary source of truth for shared run resources.
Existing CLI params remain available as fallbacks when a resource is omitted
from the YAML.

### Legacy compatibility

The parser still accepts the older list-style YAML contract as a compatibility
shim. That path is kept for low-risk transition only. The supported public
contract is the hierarchical YAML shown above.

## Example Samplesheets

Committed examples and templates:

- minimal smoke test: [`assets/samplesheet.example.yaml`](assets/samplesheet.example.yaml)
- real RNA template: [`assets/samplesheet.real_rna.template.yaml`](assets/samplesheet.real_rna.template.yaml)
- real DNA example: [`assets/samplesheet.dna.RealDATAexample.yaml`](assets/samplesheet.dna.RealDATAexample.yaml)
- real multiome example: [`assets/test_realdata/samplesheet.real_multiome.yaml`](assets/test_realdata/samplesheet.real_multiome.yaml)

## Outputs

RNA core publishes:

- `tagging/`
- `split/`
- `usam/`
- `align/`
- `pipeline_info/`

DNA core publishes:

- `dna_tagging/`
- `dna_split/`
- `dna_align/`
- `dna_dedup/`
- `dna_nodup/`
- `dna_coverage/`
- `pipeline_info/`

Optional downstream publishes:

- `shared_stage/sc_process_stage/`
- `shared_stage/sc_process_run/`

## Repo Layout

Key owned paths:

- [`main.nf`](main.nf)
- [`workflows/treseq.nf`](workflows/treseq.nf)
- [`subworkflows/local/rna_core.nf`](subworkflows/local/rna_core.nf)
- [`subworkflows/local/dna_core.nf`](subworkflows/local/dna_core.nf)
- [`subworkflows/local/optional_sc_process.nf`](subworkflows/local/optional_sc_process.nf)
- [`modules/local/`](modules/local/)
- [`scripts/core_runtime/`](scripts/core_runtime/)
- [`docs/architecture/implemented_pipeline.md`](docs/architecture/implemented_pipeline.md)

## Additional Runtime Profiles

Optional local dev micromamba/conda overlay:

```bash
micromamba env create -f envs/first_slice.yml
micromamba activate tres
nextflow run . -profile test,conda_dev --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Optional Docker overlay for the retained smoke path:

```bash
docker build -f docker/first_slice.Dockerfile -t tresflow-first-slice:py312 .
nextflow run . -profile test,docker --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Neither overlay removes the current pinned host Codon/Seq requirement.

## Validation Status

Current validated boundaries on this server:

- RNA through `ALIGN_RNA`
- DNA through `BAM_COVERAGE_DNA`
- optional shared staging as a separate downstream boundary

Current acceptance criteria:

- `-profile test` succeeds through the retained minimal smoke path
- `-profile test_real_rna` succeeds through real RNA `ALIGN_RNA`
- real DNA validation succeeds through the current real DNA boundary on this server
- optional multiome staging and optional `sc_process.py` remain real-data-only paths

## Troubleshooting

Verify pinned Codon/Seq before any run:

```bash
bin/check_codon_seq_host.sh
```

Install pinned Codon if needed:

```bash
bash scripts/install_codon_0.16.3.sh
```

Known runtime notes:

- `bamCoverage` can be long-running on this server; treat that as a runtime/performance characteristic unless there is a distinct configuration or contract error.
- The optional `sc_process.py` path is not part of the mandatory core workflow.
- The optional downstream runtime gets past SnapATAC annotation resolution using the local cache, but later downstream analysis issues remain outside the core workflow contract.
