# TrESFlow

TrESFlow is a Nextflow DSL2 wrapper around the read-only source material in `upstream/source_scripts/`. The current implementation is intentionally small: it parses a single YAML samplesheet and runs the first RNA tagging slice by wrapping the upstream `Tag.codon` and `Tag_UMI.codon` steps.

## Current slice

- `TAG_RNA_SAMPLE_BARCODE` wraps `upstream/source_scripts/Tag.codon` for the RNA sample-barcode step.
- `TAG_RNA_UMI` wraps `upstream/source_scripts/Tag_UMI.codon` for the RNA UMI step.
- `-profile test` uses lightweight mock wrappers so the pipeline runs end-to-end without Codon.
- Real runs keep the business logic in the upstream scripts and require host-installed Codon plus the Seq plugin.
- `envs/first_slice.yml` is the source of truth for software requirements around the currently implemented RNA tagging slice.
- An optional `docker` profile containerizes the current Python wrapper steps for the smoke-test path only.

## Layout

- `main.nf`
- `nextflow.config`
- `conf/base.config`
- `conf/test.config`
- `workflows/treseq.nf`
- `subworkflows/local/initial_rna_tagging.nf`
- `modules/local/tag_rna_sb/main.nf`
- `modules/local/tag_rna_umi/main.nf`
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
```

The current slice only accepts `modality: rna`. DNA and downstream alignment/splitting steps are intentionally not wired yet.

## Running the test profile

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Expected tagged outputs land under `results/test/tagging/`, and Nextflow report artifacts land under `results/test/pipeline_info/`. The trace artifact is the standard Nextflow tabular trace file, written as `execution_trace.tsv` in that directory.

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

Under `-profile docker`, those processes run in the local image `tresflow-first-slice:py312`, built from `docker/first_slice.Dockerfile`.
This does not make real non-mock execution portable, because Codon and Seq are still outside Docker in the current implementation.

## Real Mode Host Prerequisites

Real non-mock execution currently requires the following on the host:

- native Linux or macOS
- Codon CLI installed via Exaloop's installer script and available on `PATH`
- Seq plugin installed separately from the platform-specific `0.11.4` release tarball so that
  `${HOME}/.codon/lib/codon/plugins/seq/plugin.toml` exists on disk
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

Relevant entries from `envs/first_slice.yml` for the implemented RNA tagging path:

- Launcher environment when you manually `micromamba activate tres`:
  `nextflow`, `openjdk`
- Current process runtime for [`modules/local/tag_rna_sb/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_sb/main.nf) and [`modules/local/tag_rna_umi/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_umi/main.nf):
  `python` / `cpython`

Not provided by `envs/first_slice.yml`:

- `codon`
- Seq plugin files under `${HOME}/.codon/lib/codon/plugins/seq`

Current process usage:

- `TAG_RNA_SAMPLE_BARCODE` runs [`bin/run_tag.py`](/Users/aannan/GitAA/TrESFlow/bin/run_tag.py).
  Today it only requires Python's standard library in both `mock` and `real` mode.
  In `real` mode it also shells out to host-provided `codon` with `-plugin seq`, which is not currently supplied by `envs/first_slice.yml`.
- `TAG_RNA_UMI` runs [`bin/run_tag_umi.py`](/Users/aannan/GitAA/TrESFlow/bin/run_tag_umi.py).
  It has the same dependency pattern: Python stdlib plus host-provided `codon` with the Seq plugin for `real` mode.

Current Docker process coverage:

- `TAG_RNA_SAMPLE_BARCODE` is containerized for the smoke-test path
- `TAG_RNA_UMI` is containerized for the smoke-test path
- Codon and Seq are not containerized in this pass

Dependencies present in `envs/first_slice.yml` but currently unused by the implemented RNA tagging slice:

- Alignment and downstream genomics tools:
  `samtools`, `star`, `bwa-mem2`, `bedtools`, `trim-galore`, `cutadapt`, `fastqc`, `multiqc`, `deeptools`, `ucsc-bedgraphtobigwig`
- Downstream analysis stack for later steps such as `sc_process.py`:
  `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, `matplotlib`, `matplotlib-venn`, `upsetplot`, pip `snapatac2`, pip `MACS3`
- Test and developer tooling not used by the current runtime path:
  `pytest`, `pytest-timeout`
- Other packages currently unused by the implemented RNA tagging slice:
  `pyyaml`, `coreutils`, `parallel`, `pigz`, `pbzip2`

## Docker Status

Supported via Docker now:

- An optional `docker` profile exists for the current smoke-test path.
- It containerizes only the Python wrapper processes using the local image `tresflow-first-slice:py312`.
- `docker + test` is the currently supported portable smoke-test path.

What remains host-dependent right now:

- `codon` for real non-mock execution of the two wrapped upstream steps
- Seq plugin installation under `${HOME}/.codon/lib/codon/plugins/seq`
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
- Real non-mock execution is still not fully portable because the wrapped steps call `codon -plugin seq`, and Codon plus the Seq plugin are still external host prerequisites.

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
