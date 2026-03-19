# TrESFlow

TrESFlow is a Nextflow DSL2 wrapper around the read-only source material in `upstream/source_scripts/`. The current implementation is intentionally small: it parses a single YAML samplesheet and runs the current RNA-only slice by wrapping the upstream `Tag.codon`, `Tag_UMI.codon`, `Tag_Lig3.codon`, and `trim_galore` steps.

## Current slice

- `TAG_RNA_SAMPLE_BARCODE` wraps `upstream/source_scripts/Tag.codon` for the RNA sample-barcode step.
- `TAG_RNA_UMI` wraps `upstream/source_scripts/Tag_UMI.codon` for the RNA UMI step.
- `TAG_RNA_CELL_BARCODE` wraps `upstream/source_scripts/Tag_Lig3.codon` for the RNA cell-barcode step.
- `TRIM_RNA_FASTQS` wraps the immediate upstream RNA `trim_galore` step after CB tagging.
- `-profile test` uses lightweight mock wrappers so the pipeline runs end-to-end without Codon.
- Real runs keep the business logic in the upstream scripts and require host-installed Codon, the Seq plugin, and `trim_galore`.
- `envs/first_slice.yml` is the source of truth for software requirements around the currently implemented RNA-only slice.
- An optional `docker` profile containerizes the current Python wrapper steps for the smoke-test path only.
- The next upstream RNA step after `Tag_Lig3` is `trim_galore`. `Split_ReadsV2`, `FqToSAM`, and `AlignRNA.sh` remain intentionally out of scope for this slice.

## Layout

- `main.nf`
- `nextflow.config`
- `conf/base.config`
- `conf/test.config`
- `workflows/treseq.nf`
- `subworkflows/local/initial_rna_tagging.nf`
- `modules/local/tag_rna_sb/main.nf`
- `modules/local/tag_rna_umi/main.nf`
- `modules/local/tag_rna_cell_barcode/main.nf`
- `modules/local/trim_rna_fastqs/main.nf`
- `assets/samplesheet.example.yaml`
- `assets/testdata/`

## Required CLI params

- `--samplesheet`
- `--outdir`

## Optional CLI params

- `--upstream_dir`
  Default: `./upstream/source_scripts`

## Samplesheet schema

```yaml
samples:
  - id: sample_id
    modality: rna
    reads:
      i1: path/to/I1.fastq.gz
      r1: path/to/R1.fastq.gz
      r2: path/to/R2.fastq.gz
    barcodes:
      sample:
        whitelist: path/to/sample_whitelist.txt
        bc_len: 4
        bc_start: 0
        hd: 1
        tag: SB
        first_pass: first_pass
        reverse_complement: true
      umi:
        bc_len: 10
        bc_start: 4
        tag: UM
      cell:
        whitelist: path/to/cell_whitelist.txt
        bc_len: 8
        hd: 1
        tag: CB
```

The current slice only accepts `modality: rna`. DNA and downstream alignment/splitting steps are intentionally not wired yet.

## Running the test profile

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Expected tagging and trimming outputs land under `results/test/tagging/`, and Nextflow report artifacts land under `results/test/pipeline_info/`. The trace artifact is the standard Nextflow tabular trace file, written as `execution_trace.tsv` in that directory.
The upstream launcher deletes the untrimmed CB FASTQs after `trim_galore`; this slice keeps them published as intermediates and advances on the trimmed `_val_1` / `_val_2` outputs.

## Current RNA Step Map

The upstream RNA order currently relevant to this repo is:

1. `Tag.codon` for sample barcode (`SB`)
2. `Tag_UMI.codon` for UMI (`UM`)
3. `Tag_Lig3.codon` for cell barcode (`CB`) and `RG`
4. `trim_galore`
5. `Split_ReadsV2.codon` in `rna` mode
6. `FqToSAM.codon`
7. `AlignRNA.sh`

This repo now implements steps 1 through 4 only. Under normal execution they are real by default. Under `-profile test`, all four wrapped RNA steps use mock behavior.

## Micromamba and Conda

If you want the launcher environment to match the checked-in dependency spec, create and activate the micromamba environment defined in `envs/first_slice.yml`:

```bash
micromamba env create -f envs/first_slice.yml
micromamba activate tres
```

That environment includes `nextflow` and `openjdk`, so it is a convenient way to launch the current code path locally.
It does not install Codon or the Seq plugin.

An optional Nextflow profile is also available for task-level environment management using Nextflow's official `conda` support with micromamba enabled:

```bash
nextflow run . -profile test,conda_dev --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

`conda_dev` is for local development only. It is not enabled by default, and it does not replace Docker as the preferred portable execution target.
Today it fully covers the mock `-profile test` path. Real mode still depends on host-installed Codon and the Seq plugin.

## What Works Today

### Host + Test

Supported on native Linux and macOS hosts:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

This path uses the bundled mock wrappers. It does not require Codon or Seq.

### Micromamba Dev + Test

Supported on native Linux and macOS hosts:

```bash
micromamba env create -f envs/first_slice.yml
micromamba activate tres
nextflow run . -profile test,conda_dev --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

This path uses Nextflow's official `conda` support with `conda.useMicromamba = true`. It is useful for local development and keeps the launcher and task Python environment aligned with `envs/first_slice.yml`.

### Docker + Test

Supported on native Linux and macOS hosts with Docker available:

Build the small first-slice image once:

```bash
docker build -f docker/first_slice.Dockerfile -t tresflow-first-slice:py312 .
```

Then run:

```bash
nextflow run . -profile test,docker --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

This is the portable smoke-test path for the currently implemented first slice.
Only the current Python wrapper processes are containerized in this pass:

- [`TAG_RNA_SAMPLE_BARCODE`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_sb/main.nf)
- [`TAG_RNA_UMI`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_umi/main.nf)
- [`TAG_RNA_CELL_BARCODE`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_cell_barcode/main.nf)
- [`TRIM_RNA_FASTQS`](/Users/aannan/GitAA/TrESFlow/modules/local/trim_rna_fastqs/main.nf)

Under `-profile docker`, those processes run in the local image `tresflow-first-slice:py312`, built from `docker/first_slice.Dockerfile`.
This does not make real non-mock execution portable, because Codon and Seq are still outside Docker in the current implementation.

## Real Mode Host Prerequisites

Real non-mock execution currently requires the following on the host:

- native Linux or macOS
- Codon CLI installed via Exaloop's installer script and available on `PATH`
- Seq plugin installed separately from the platform-specific `0.11.4` release tarball so that
  `${HOME}/.codon/lib/codon/plugins/seq/plugin.toml` exists on disk
- `trim_galore` available on `PATH` for the RNA trimming step
- the read-only upstream scripts under `upstream/source_scripts/`

Pinned Seq example for current documentation:

- Seq plugin version: `0.11.4`
- plugin metadata on this host reports compatibility: `supported = ">=0.18.2"`
- Codon installed on this host reports version: `0.19.6`

Example verification commands:

```bash
codon --version
sed -n '1,20p' "${HOME}/.codon/lib/codon/plugins/seq/plugin.toml"
bin/check_codon_seq_host.sh
```

Example real-mode run:

```bash
nextflow run . --samplesheet /path/to/samplesheet.yaml --outdir /path/to/results
```

For real execution, `envs/first_slice.yml` is not enough by itself. Codon and Seq remain documented host prerequisites until they are containerized.

## Dependency Mapping for the Current Slice

Relevant entries from `envs/first_slice.yml` for the implemented RNA-only path:

- Launcher environment when you manually `micromamba activate tres`:
  `nextflow`, `openjdk`
- Current wrapper runtime for [`modules/local/tag_rna_sb/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_sb/main.nf), [`modules/local/tag_rna_umi/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_umi/main.nf), [`modules/local/tag_rna_cell_barcode/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_cell_barcode/main.nf), and [`modules/local/trim_rna_fastqs/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/trim_rna_fastqs/main.nf):
  `python` / `cpython`
- Current real trimming runtime for [`modules/local/trim_rna_fastqs/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/trim_rna_fastqs/main.nf):
  `trim-galore`, `cutadapt`

Not provided by `envs/first_slice.yml`:

- `codon`
- Seq plugin files under `${HOME}/.codon/lib/codon/plugins/seq`

Current process usage:

- `TAG_RNA_SAMPLE_BARCODE` runs [`bin/run_tag.py`](/Users/aannan/GitAA/TrESFlow/bin/run_tag.py).
  Today it only requires Python's standard library in both `mock` and `real` mode.
  In `real` mode it also shells out to host-provided `codon` with `-plugin seq`, which is not currently supplied by `envs/first_slice.yml`.
- `TAG_RNA_UMI` runs [`bin/run_tag_umi.py`](/Users/aannan/GitAA/TrESFlow/bin/run_tag_umi.py).
  It has the same dependency pattern: Python stdlib plus host-provided `codon` with the Seq plugin for `real` mode.
- `TAG_RNA_CELL_BARCODE` runs [`bin/run_tag_lig3.py`](/Users/aannan/GitAA/TrESFlow/bin/run_tag_lig3.py).
  It also uses Python stdlib in `mock` mode and shells out to host-provided `codon` with the Seq plugin in `real` mode.
- `TRIM_RNA_FASTQS` runs [`bin/run_trim_galore.py`](/Users/aannan/GitAA/TrESFlow/bin/run_trim_galore.py).
  In `mock` mode it gzip-copies the CB-tagged FASTQs to the expected trim_galore `_val_1` / `_val_2` outputs.
  In `real` mode it shells out to `trim_galore` with the upstream launcher settings `--quality 10 --gzip --length 20 --paired`.

Current Docker process coverage:

- `TAG_RNA_SAMPLE_BARCODE` is containerized for the smoke-test path
- `TAG_RNA_UMI` is containerized for the smoke-test path
- `TAG_RNA_CELL_BARCODE` is containerized for the smoke-test path
- `TRIM_RNA_FASTQS` is containerized for the smoke-test path
- Codon and Seq are not containerized in this pass

Dependencies present in `envs/first_slice.yml` but currently unused by the implemented RNA-only slice:

- Alignment and downstream genomics tools:
  `samtools`, `star`, `bwa-mem2`, `bedtools`, `fastqc`, `multiqc`, `deeptools`, `ucsc-bedgraphtobigwig`
- Downstream analysis stack for later steps such as `sc_process.py`:
  `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, `matplotlib`, `matplotlib-venn`, `upsetplot`, pip `snapatac2`, pip `MACS3`
- Test and developer tooling not used by the current runtime path:
  `pytest`, `pytest-timeout`
- Other packages currently unused by the implemented RNA-only slice:
  `pyyaml`, `coreutils`, `parallel`, `pigz`, `pbzip2`

## Docker Status

Supported via Docker now:

- An optional `docker` profile exists for the current smoke-test path.
- It containerizes only the Python wrapper processes using the local image `tresflow-first-slice:py312`.
- `docker + test` is the currently supported portable smoke-test path.

What remains host-dependent right now:

- `codon` for real non-mock execution of the three wrapped upstream Codon RNA steps
- Seq plugin installation under `${HOME}/.codon/lib/codon/plugins/seq`
- `trim_galore` for real non-mock execution of the RNA trim step unless you choose the optional `-profile conda_dev`
- The read-only upstream scripts under `upstream/source_scripts/`
- The local launcher environment unless you choose either:
  `micromamba activate tres`
  or the optional `-profile conda_dev`

Containerization note for the upcoming Docker work:

- `envs/first_slice.yml` is now a portable dependency manifest and should be treated as the package baseline for the first-slice image.
- Docker remains the preferred long-term portability target for Linux and macOS users running containers.
- A Docker image that fully supports real execution of the current slice will still need to account for both Codon and the Seq plugin, because those requirements are outside the checked-in env file today.

Current limitation for real mode:

- `-profile docker` only makes the mock smoke test portable today.
- Real non-mock execution is still not fully portable because the wrapped steps call `codon -plugin seq`, Codon plus the Seq plugin are still external host prerequisites, and the Docker smoke-test image does not include `trim_galore`.

## Acceptance Criteria

The current RNA-only slice is accepted when:

- `nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test` completes successfully
- `nextflow run . -profile test,docker --samplesheet assets/samplesheet.example.yaml --outdir results/test` completes successfully
- `results/test/tagging/` contains:
  `test_rna.sample_barcode.R1.fastq`
  `test_rna.sample_barcode.R2.fastq`
  `test_rna.sample_barcode.counts.tsv`
  `test_rna.sample_barcode.stats.tsv`
  `test_rna.sample_barcode_umi.R1.fastq`
  `test_rna.sample_barcode_umi.R2.fastq`
  `test_rna.umi.counts.tsv`
  `test_rna.sample_barcode_umi_cell.R1.fastq`
  `test_rna.sample_barcode_umi_cell.R2.fastq`
  `test_rna.cell.counts.tsv`
  `test_rna.tag_records.tsv`
  `test_rna.cell.stats_L1.tsv`
  `test_rna.cell.stats_L2.tsv`
  `test_rna.cell.stats_L3.tsv`
  `test_rna.sample_barcode_umi_cell.R1_val_1.fq.gz`
  `test_rna.sample_barcode_umi_cell.R2_val_2.fq.gz`
- `results/test/pipeline_info/` contains the configured Nextflow report artifacts

Exact test command:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

## Troubleshooting

To verify the host prerequisites for real mode:

```bash
bin/check_codon_seq_host.sh
```

That check verifies:

- `codon` is on `PATH`
- Codon reports a version
- Seq plugin metadata exists where Codon expects it by default
- Seq plugin metadata exposes a version and supported Codon range
