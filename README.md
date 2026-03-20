# TrESFlow

TrESFlow is a Nextflow DSL2 pipeline that owns the **implemented core runtime**
inside this repo.

- Core workflow:
  - RNA through `ALIGN_RNA`
  - DNA through `BAM_COVERAGE_DNA`
- Optional downstream:
  - shared staging
  - one optional `sc_process.py` call

Architecture/DAG view:
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

`RUN_SC_PROCESS` is explicit and optional. The core pipeline is complete without
it, and the workflow design allows only one `sc_process.py` call.

## Quick Start

Mock RNA core:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Mock DNA core:

```bash
nextflow run . -profile test_dna --samplesheet assets/samplesheet.dna.example.yaml --outdir results/test_dna
```

Mock multiome plus shared staging:

```bash
nextflow run . -profile test_multiome --samplesheet assets/samplesheet.multiome.example.yaml --outdir results/test_multiome
```

Real RNA validation:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet assets/test_realdata/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna \
  --rna_ref_base_dir /path/to/reference_base \
  --rna_align_species human
```

Real DNA validation:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.dna.RealDATAexample.yaml \
  --outdir results/test_dna_real \
  --dna_bwa_reference /path/to/bwa_index_prefix \
  --dna_blacklist_bed /path/to/blacklist.bed \
  --dna_effective_genome_size 2913022398
```

Optional shared staging only:

```bash
nextflow run . \
  -profile test_multiome \
  --samplesheet assets/samplesheet.multiome.example.yaml \
  --outdir results/test_multiome_stage \
  --stage_sc_process_inputs true
```

Optional downstream analysis:

```bash
nextflow run . \
  --samplesheet assets/test_realdata/samplesheet.real_multiome.yaml \
  --outdir results/test_real_multiome_sc_process \
  --rna_ref_base_dir /path/to/reference_base \
  --rna_align_species human \
  --dna_bwa_reference /path/to/bwa_index_prefix \
  --dna_blacklist_bed /path/to/blacklist.bed \
  --dna_effective_genome_size 2913022398 \
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
The read-only [`upstream/source_scripts/`](upstream/source_scripts/) tree remains
for provenance and for the optional `sc_process.py` entrypoint only.

## Inputs and Params

### Required for every run

- `--samplesheet`
- `--outdir`

Pipeline input is a single YAML samplesheet under `params.samplesheet`.

### Shared samplesheet contract

Top-level keys:

- `library_name`
- `sb_group_map`
- `dna_mo_map` when any DNA samples are present

The parser accepts `modality: rna` and `modality: dna` rows in the same file.

`sb_group_map` is shared by RNA and DNA. It is:

- the sample-barcode grouping TSV used by `Split_ReadsV2.codon`
- the source of truth for the effective DNA sample-barcode whitelist

The current real-data examples still use the legacy on-disk filename
`sb_map_RNA.tsv`, but the pipeline treats it as the shared `sb_group_map`.

### Required for real RNA runs

- `--rna_ref_base_dir`
- `--rna_align_species`

Reference layout:

- human: `<ref_base_dir>/GRCh38_TrES/star` and `<ref_base_dir>/hg38.chrom.sizes`
- mouse: `<ref_base_dir>/GRCm39_TrES/star` and `<ref_base_dir>/mm39.chrom.sizes`

### Required for real DNA runs

- `--dna_bwa_reference`
- `--dna_blacklist_bed`
- `--dna_effective_genome_size`

`--dna_bwa_reference` is a bwa-mem2 index prefix. The following sidecars must
exist: `.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`.

### Optional downstream controls

- `--stage_sc_process_inputs`
- `--run_sc_process`
- `--runtime_snap_data_dir`
- `--optional_sc_process_script`

`RUN_SC_PROCESS` binds `SNAP_DATA_DIR` explicitly. By default on this server it
uses `/home/annan/.cache/snapatac2`, so `snap.genome.hg38` and
`snap.genome.mm39` resolve from the local SnapATAC cache instead of requiring a
new user-supplied annotation file by default.

### Runtime and layout overrides

- `--core_scripts_dir`
- `--runtime_bin_dir`
- individual `--runtime_*` params
- `--runtime_codon`
- `--codon_home`
- `--max_cpus`

`--max_cpus` defaults to `40`.

## Example Samplesheets

Committed examples:

- RNA mock: [`assets/samplesheet.example.yaml`](assets/samplesheet.example.yaml)
- DNA mock: [`assets/samplesheet.dna.example.yaml`](assets/samplesheet.dna.example.yaml)
- Multiome mock: [`assets/samplesheet.multiome.example.yaml`](assets/samplesheet.multiome.example.yaml)
- Real RNA template: [`assets/samplesheet.real_rna.template.yaml`](assets/samplesheet.real_rna.template.yaml)
- Real multiome example: [`assets/test_realdata/samplesheet.real_multiome.yaml`](assets/test_realdata/samplesheet.real_multiome.yaml)

Committed test assets live under:

- [`assets/testdata/`](assets/testdata/)
- [`assets/test_realdata/`](assets/test_realdata/)

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

## Additional Profiles

Local dev micromamba/conda support:

```bash
micromamba env create -f envs/first_slice.yml
micromamba activate tres
nextflow run . -profile test,conda_dev --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Docker smoke path:

```bash
docker build -f docker/first_slice.Dockerfile -t tresflow-first-slice:py312 .
nextflow run . -profile test,docker --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Neither of those removes the current pinned host Codon/Seq requirement.

## Validation Status

Current validated boundaries on this server:

- RNA through `ALIGN_RNA`
- DNA through `BAM_COVERAGE_DNA`
- optional shared staging as a separate downstream boundary

Current acceptance criteria:

- `-profile test` succeeds through mocked RNA `ALIGN_RNA`
- `-profile test_dna` succeeds through mocked DNA `BAM_COVERAGE_DNA`
- `-profile test_multiome` succeeds through shared staging
- `-profile test_real_rna` succeeds through real RNA `ALIGN_RNA`
- real DNA validation succeeds through direct `_NoDup.bam` outputs and reaches the coverage wrapper with the established server runtime contract

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

- `bamCoverage` can be long-running on this server; treat it as a runtime/performance characteristic unless there is a distinct configuration or contract error.
- The optional `sc_process.py` path is not part of the mandatory core workflow.
- The optional downstream runtime now gets past SnapATAC annotation resolution using the local cache, but later downstream analysis issues remain outside the core workflow contract.
