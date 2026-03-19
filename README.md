# TrESFlow

TrESFlow is a Nextflow DSL2 wrapper around the read-only source material in `upstream/source_scripts/`. The current implementation is intentionally small: it parses a single YAML samplesheet and runs the current RNA-only slice by wrapping the upstream `Tag.codon`, `Tag_UMI.codon`, `Tag_Lig3.codon`, `trim_galore`, `Split_ReadsV2.codon` RNA-mode, `FqToSAM.codon`, and optionally `AlignRNA.sh`.

## Current slice

- `TAG_RNA_SAMPLE_BARCODE` wraps `upstream/source_scripts/Tag.codon` for the RNA sample-barcode step.
- `TAG_RNA_UMI` wraps `upstream/source_scripts/Tag_UMI.codon` for the RNA UMI step.
- `TAG_RNA_CELL_BARCODE` wraps `upstream/source_scripts/Tag_Lig3.codon` for the RNA cell-barcode step.
- `TRIM_RNA_FASTQS` wraps the immediate upstream RNA `trim_galore` step after CB tagging.
- `SPLIT_RNA_READS` wraps the immediate upstream RNA `Split_ReadsV2.codon` step in `rna` mode after trimming.
- `FQ_TO_SAM` wraps the immediate upstream RNA `FqToSAM.codon` step after read splitting.
- `ALIGN_RNA` wraps the immediate upstream RNA `AlignRNA.sh` step after grouped `.usam` generation when `--run_align_rna true`.
- `-profile test` uses lightweight mock wrappers for the wrapped RNA steps, but the pipeline still enforces host Codon `0.16.3` and Seq `0.11.3` before any run starts.
- Every pipeline run requires host-installed Codon `0.16.3` and Seq `0.11.3`. Real preprocessing runs also require `trim_galore`. Real RNA alignment also requires `STAR`, `samtools`, and `bedGraphToBigWig`.
- `envs/first_slice.yml` is the source of truth for software requirements around the currently implemented RNA-only slice.
- A separate `test_real_rna` profile is available for external/local validation data through the current preprocessing boundary, with optional extension through `AlignRNA.sh`.
- An optional `docker` profile containerizes the current Python wrapper steps for the smoke-test path only.
- Pinning the real-mode toolchain does not change the mock `-profile test` or `-profile test,docker` paths.
- The next upstream RNA step after `FqToSAM.codon` is `AlignRNA.sh`. Downstream `samtools view` to `.ubam` plus `sc_process.py` remain intentionally out of scope for this slice.

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
- `modules/local/split_rna_reads/main.nf`
- `modules/local/fq_to_sam/main.nf`
- `modules/local/align_rna/main.nf`
- `assets/samplesheet.example.yaml`
- `assets/samplesheet.real_rna.template.yaml`
- `assets/testdata/`
- `scripts/install_codon_0.16.3.sh`

## Required CLI params

- `--samplesheet`
- `--outdir`

## Optional CLI params

- `--upstream_dir`
  Default: `./upstream/source_scripts`
- `--run_align_rna`
  Default: `false`
- `--rna_ref_base_dir`
  Required when `--run_align_rna true`. Must contain either `GRCh38_TrES/star` plus `hg38.chrom.sizes`, or `GRCm39_TrES/star` plus `mm39.chrom.sizes`, depending on species.
- `--rna_align_species`
  Required when `--run_align_rna true`. Supported values: `human`, `mouse`

## Samplesheet schema

```yaml
library_name: run_library_name
rna_sb_group_map: path/to/rna_sb_group_map.tsv

samples:
  - id: sample_id
    modality: rna
    reads:
      i1: path/to/I1.fastq.gz
      r1: path/to/R1.fastq.gz
      r2: path/to/R2.fastq.gz
    barcodes:
      sample:
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

The current slice only accepts `modality: rna`. DNA plus downstream alignment steps are intentionally not wired yet.
The RNA SB-group map is the launcher-style TSV used by `Split_ReadsV2.codon` in `rna` mode: `sample<TAB>sb_group<TAB>sb_bc`.
For the current RNA workflow it is also the single source of truth for the sample-barcode whitelist passed into `Tag.codon`.
Blank lines, `#` comments, and a literal `sample sb_group sb_bc` header row are ignored.
The pipeline fails if the map has no rows for a sample or if the same `sb_bc` is assigned to multiple groups for that sample.

## Running the test profile

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Expected tagging and trimming outputs land under `results/test/tagging/`, split outputs land under `results/test/split/`, grouped unmapped SAM outputs land under `results/test/usam/`, optional mock alignment outputs land under `results/test/align/` when `--run_align_rna true`, and Nextflow report artifacts land under `results/test/pipeline_info/`. The trace artifact is the standard Nextflow tabular trace file, written as `execution_trace.tsv` in that directory.
The upstream launcher deletes the untrimmed CB FASTQs after `trim_galore`; this slice keeps them published as intermediates and advances on the trimmed `_val_1` / `_val_2` outputs.
`Split_ReadsV2.codon` has one ambiguity in its comments versus examples: the code comments discuss dropping an injected leading base from `SB`, but the upstream RNA map example uses full `SB` strings. This repo follows the script's actual lookup behavior: raw `SB` match first, then drop-first fallback.

## Real RNA Validation

The real-data validation path is intentionally external/local. No real RNA fixture is committed to this repo.

Minimum real-input contract for the current RNA-only workflow:

- top-level `library_name`
- top-level `rna_sb_group_map`
- per-sample `id` matching the `sample` column in `rna_sb_group_map` for the groups you want emitted.
  For the current real-mode example, that sample id is `day15`.
- per-sample `modality: rna`
- per-sample `reads.i1`
- per-sample `reads.r1`
- per-sample `reads.r2`
- per-sample `barcodes.sample.bc_len`
- per-sample `barcodes.sample.bc_start`
- per-sample `barcodes.sample.hd`
- per-sample `barcodes.sample.tag` or default `SB`
- per-sample `barcodes.sample.first_pass` or default `first_pass`
- per-sample `barcodes.sample.reverse_complement` or default `true`
- per-sample `barcodes.umi.bc_len`
- per-sample `barcodes.umi.bc_start`
- per-sample `barcodes.umi.tag` or default `UM`
- per-sample `barcodes.cell.whitelist`
- per-sample `barcodes.cell.bc_len`
- per-sample `barcodes.cell.hd`
- per-sample `barcodes.cell.tag` or default `CB`

Recommended external/local layout:

```text
/path/to/real_rna_validation/
  samplesheet.real_rna.yaml
  sb_map_RNA.tsv
  ligation_barcode_whitelist.txt
  day15_I1.fq.gz
  day15_R1.fq.gz
  day15_R2.fq.gz
```

Use [`samplesheet.real_rna.template.yaml`](/Users/aannan/GitAA/TrESFlow/assets/samplesheet.real_rna.template.yaml) as the committed template. The real files above are expected to remain external/local.
The template keeps `samples[0].id: day15` because the sample id must match the first `sample` column in `rna_sb_group_map`.

Exact run command for host-based real-RNA preprocessing validation:

```bash
nextflow run . -profile test_real_rna --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml --outdir results/test_real_rna
```

Exact run command for host-based real-RNA alignment validation through `AlignRNA.sh`:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna_align \
  --run_align_rna true \
  --rna_ref_base_dir /mnt/dataFast/ahrmad \
  --rna_align_species human
```

Optional variant if you want Nextflow to provision `trim_galore` and Python task dependencies through the checked-in conda environment while still using the globally required host Codon `0.16.3` plus Seq `0.11.3`:

```bash
nextflow run . -profile test_real_rna,conda_dev --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml --outdir results/test_real_rna
```

Expected outputs through the current preprocessing boundary:

- `outdir/tagging/` with SB, UM, CB, and trim outputs
- `outdir/split/` with one FASTQ pair and one `SAM_RG_Header_*.tsv` per RNA SB group in `rna_sb_group_map`
- `outdir/usam/` with one `<sample>_<group>_tagged.usam` per RNA SB group in `rna_sb_group_map`
- `outdir/pipeline_info/` with the configured Nextflow report artifacts

Expected additional outputs when `--run_align_rna true`:

- `outdir/align/` with one `<sample>_<group>.Solo.outGeneFull/` STARsolo directory per grouped `.usam`
- `outdir/align/` with one `<sample>_<group>.filtered_cells.bam` per grouped `.usam`
- `outdir/align/` with stranded and unstranded bigWig tracks per grouped `.usam` when `AlignRNA.sh` emits them

## Current RNA Step Map

The upstream RNA order currently relevant to this repo is:

1. `Tag.codon` for sample barcode (`SB`)
2. `Tag_UMI.codon` for UMI (`UM`)
3. `Tag_Lig3.codon` for cell barcode (`CB`) and `RG`
4. `trim_galore`
5. `Split_ReadsV2.codon` in `rna` mode
6. `FqToSAM.codon`
7. `AlignRNA.sh`

This repo now implements steps 1 through 7. Steps 1 through 6 remain the validated preprocessing boundary. Step 7 is opt-in via `--run_align_rna true`. Under `-profile test`, all seven wrapped RNA steps use mock behavior.

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
Today it fully covers the mock `-profile test` task runtime. Every pipeline run still depends on host-installed Codon `0.16.3` and Seq `0.11.3`.

## What Works Today

### Host + Test

Supported on native Linux and macOS hosts:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

This path uses the bundled mock wrappers for the RNA steps, but startup still fails unless host Codon `0.16.3` and Seq `0.11.3` are installed.

### Host + Real RNA Validation

Supported on native Linux and macOS hosts with external/local RNA inputs:

```bash
nextflow run . -profile test_real_rna --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml --outdir results/test_real_rna
```

This path runs the current RNA-only workflow in real mode through grouped `FQ_TO_SAM` outputs.
It requires host-installed Codon `0.16.3`, Seq `0.11.3`, and `trim_galore`.
It expects the real input files to remain external/local and not committed to this repo.

### Host + Real RNA Alignment

Supported on this server's host environment when the upstream STAR reference base exists and you supply `--rna_ref_base_dir` plus `--rna_align_species`:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna_align \
  --run_align_rna true \
  --rna_ref_base_dir /mnt/dataFast/ahrmad \
  --rna_align_species human
```

This path runs the current RNA-only workflow through grouped `AlignRNA.sh` outputs.
It requires the same preprocessing prerequisites plus host `STAR`, `samtools`, and `bedGraphToBigWig`.
It also requires the reference-base contract expected by `AlignRNA.sh`:
`<ref_base_dir>/GRCh38_TrES/star` with `<ref_base_dir>/hg38.chrom.sizes` for `human`, or `<ref_base_dir>/GRCm39_TrES/star` with `<ref_base_dir>/mm39.chrom.sizes` for `mouse`.

### Micromamba Dev + Test

Supported on native Linux and macOS hosts:

```bash
micromamba env create -f envs/first_slice.yml
micromamba activate tres
nextflow run . -profile test,conda_dev --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

This path uses Nextflow's official `conda` support with `conda.useMicromamba = true`. It is useful for local development and keeps the launcher and task Python environment aligned with `envs/first_slice.yml`.
It still requires host Codon `0.16.3` and Seq `0.11.3` before the pipeline will start.

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

This keeps the currently implemented RNA-only slice containerized for the wrapper tasks, but it is no longer a standalone portable smoke-test path because startup still enforces host Codon `0.16.3` and Seq `0.11.3`.
Only the current Python wrapper processes are containerized in this pass:

- [`TAG_RNA_SAMPLE_BARCODE`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_sb/main.nf)
- [`TAG_RNA_UMI`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_umi/main.nf)
- [`TAG_RNA_CELL_BARCODE`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_cell_barcode/main.nf)
- [`TRIM_RNA_FASTQS`](/Users/aannan/GitAA/TrESFlow/modules/local/trim_rna_fastqs/main.nf)
- [`SPLIT_RNA_READS`](/Users/aannan/GitAA/TrESFlow/modules/local/split_rna_reads/main.nf)
- [`FQ_TO_SAM`](/Users/aannan/GitAA/TrESFlow/modules/local/fq_to_sam/main.nf)
- [`ALIGN_RNA`](/Users/aannan/GitAA/TrESFlow/modules/local/align_rna/main.nf)

Under `-profile docker`, those processes run in the local image `tresflow-first-slice:py312`, built from `docker/first_slice.Dockerfile`.
This does not make any execution mode fully portable, because Codon `0.16.3` and Seq `0.11.3` are still required on the host in the current implementation.

## Real Mode Host Prerequisites

Every pipeline run currently requires the following on the host:

- native Linux or macOS
- Codon `0.16.3` installed under `${HOME}/.codon` and available on `PATH`
- Seq `0.11.3` installed separately under `${HOME}/.codon/lib/codon/plugins/seq`
- `trim_galore` available on `PATH` for real RNA trimming runs
- `STAR`, `samtools`, and `bedGraphToBigWig` available on `PATH` for real RNA alignment runs when `--run_align_rna true`
- the read-only upstream scripts under `upstream/source_scripts/`

For the `test_real_rna` profile specifically, you also need an external/local samplesheet plus the referenced `I1`, `R1`, `R2`, cell whitelist, and RNA SB-group map files. The real-RNA sample id must match the `sample` column in that map. For the current example, that value is `day15`.

If you already have a newer `${HOME}/.codon` install, back it up before downgrading:

```bash
mv "${HOME}/.codon" "${HOME}/.codon.backup-$(date +%Y%m%d%H%M%S)"
```

Install pinned Codon `0.16.3` with the repo-local helper:

```bash
bash scripts/install_codon_0.16.3.sh
```

Install pinned Seq `0.11.3` separately from the OS/ARCH-specific release tarball. Examples:

macOS arm64:

```bash
mkdir -p "${HOME}/.codon/lib/codon/plugins"
curl -LO https://github.com/exaloop/seq/releases/download/v0.11.3/seq-darwin-arm64.tar.gz
tar zxvf seq-darwin-arm64.tar.gz -C "${HOME}/.codon/lib/codon/plugins"
```

Linux x86_64:

```bash
mkdir -p "${HOME}/.codon/lib/codon/plugins"
curl -LO https://github.com/exaloop/seq/releases/download/v0.11.3/seq-linux-x86_64.tar.gz
tar zxvf seq-linux-x86_64.tar.gz -C "${HOME}/.codon/lib/codon/plugins"
```

Use the matching Seq `0.11.3` tarball name for other supported Linux/macOS OS/ARCH combinations.

Exact preflight command used for every pipeline run:

```bash
bin/check_codon_seq_host.sh
```

Example version checks:

```bash
codon --version
sed -n '1,20p' "${HOME}/.codon/lib/codon/plugins/seq/plugin.toml"
```

Example real-mode run:

```bash
nextflow run . --samplesheet /path/to/samplesheet.yaml --outdir /path/to/results
```

For every execution mode, `envs/first_slice.yml` is not enough by itself. It does not install Codon `0.16.3` or Seq `0.11.3`, which remain mandatory host prerequisites until they are containerized.

## Dependency Mapping for the Current Slice

Relevant entries from `envs/first_slice.yml` for the implemented RNA-only path:

- Launcher environment when you manually `micromamba activate tres`:
  `nextflow`, `openjdk`
- Current wrapper runtime for [`modules/local/tag_rna_sb/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_sb/main.nf), [`modules/local/tag_rna_umi/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_umi/main.nf), [`modules/local/tag_rna_cell_barcode/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/tag_rna_cell_barcode/main.nf), [`modules/local/trim_rna_fastqs/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/trim_rna_fastqs/main.nf), [`modules/local/split_rna_reads/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/split_rna_reads/main.nf), and [`modules/local/fq_to_sam/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/fq_to_sam/main.nf):
  `python` / `cpython`
- Current real trimming runtime for [`modules/local/trim_rna_fastqs/main.nf`](/Users/aannan/GitAA/TrESFlow/modules/local/trim_rna_fastqs/main.nf):
  `trim-galore`, `cutadapt`

Not provided by `envs/first_slice.yml`:

- Codon `0.16.3`
- Seq `0.11.3` plugin files under `${HOME}/.codon/lib/codon/plugins/seq`

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
- `SPLIT_RNA_READS` runs [`bin/run_split_reads_rna.py`](/Users/aannan/GitAA/TrESFlow/bin/run_split_reads_rna.py).
  In `mock` mode it reproduces the upstream RNA grouping behavior from the trimmed FASTQs and writes launcher-style `sample_group_R1/R2.fq.gz` plus `SAM_RG_Header_sample_group.tsv` files.
  In `real` mode it shells out to `codon -plugin seq` for `Split_ReadsV2.codon` in `rna` mode with the launcher-style RNA SB-group map.
- `FQ_TO_SAM` runs [`bin/run_fq_to_sam.py`](/Users/aannan/GitAA/TrESFlow/bin/run_fq_to_sam.py).
  In `mock` mode it reproduces the checked-in `FqToSAM.codon` behavior from split FASTQs and writes grouped unmapped SAM files.
  In `real` mode it shells out to `codon -plugin seq` for `FqToSAM.codon`, which accepts `.fq.gz` inputs directly.
- `ALIGN_RNA` runs [`upstream/source_scripts/AlignRNA.sh`](/Users/aannan/GitAA/TrESFlow/upstream/source_scripts/AlignRNA.sh) directly.
  In `mock` mode it writes placeholder STARsolo directories, filtered BAMs, and bigWig files for smoke-test coverage.
  In `real` mode it shells out to host `STAR`, `samtools`, and `bedGraphToBigWig`, using the grouped `.usam` plus the reference-base and species contract expected by the upstream shell script.

Current Docker process coverage:

- `TAG_RNA_SAMPLE_BARCODE` is containerized for the smoke-test path
- `TAG_RNA_UMI` is containerized for the smoke-test path
- `TAG_RNA_CELL_BARCODE` is containerized for the smoke-test path
- `TRIM_RNA_FASTQS` is containerized for the smoke-test path
- `SPLIT_RNA_READS` is containerized for the smoke-test path
- `FQ_TO_SAM` is containerized for the smoke-test path
- `ALIGN_RNA` is containerized for the smoke-test path
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
- `docker + test` still containerizes the current wrapper tasks, but it now also requires host Codon `0.16.3` and Seq `0.11.3` because startup enforces the pinned host toolchain globally.

What remains host-dependent right now:

- Codon `0.16.3` for real non-mock execution of the five wrapped upstream Codon RNA steps
- Seq `0.11.3` installation under `${HOME}/.codon/lib/codon/plugins/seq`
- `trim_galore` for real non-mock execution of the RNA trim step unless you choose the optional `-profile conda_dev`
- The read-only upstream scripts under `upstream/source_scripts/`
- The local launcher environment unless you choose either:
  `micromamba activate tres`
  or the optional `-profile conda_dev`

Containerization note for the upcoming Docker work:

- `envs/first_slice.yml` is now a portable dependency manifest and should be treated as the package baseline for the first-slice image.
- Docker remains the preferred long-term portability target for Linux and macOS users running containers.
- A Docker image that fully supports real execution of the current slice will still need to account for both Codon `0.16.3` and Seq `0.11.3`, because those requirements are outside the checked-in env file today.

Current limitation:

- No profile bypasses the pinned host Codon/Seq requirement.
- `-profile docker` only containerizes the wrapped task runtime. It does not remove the host requirement for Codon `0.16.3` plus Seq `0.11.3`, and the Docker smoke-test image still does not include `trim_galore`.
- Real `--run_align_rna true` execution is host-only in the current implementation, because the checked-in Docker image does not include `STAR`, `samtools`, or `bedGraphToBigWig`.

## Acceptance Criteria

The current RNA-only slice is accepted when:

- `nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test` completes successfully
- `nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test_align_mock --run_align_rna true --rna_ref_base_dir /mnt/dataFast/ahrmad --rna_align_species human` completes successfully
- `nextflow run . -profile test,docker --samplesheet assets/samplesheet.example.yaml --outdir results/test` completes successfully
- `nextflow run . -profile test,docker --samplesheet assets/samplesheet.example.yaml --outdir results/test_align_mock_docker --run_align_rna true --rna_ref_base_dir /mnt/dataFast/ahrmad --rna_align_species human` completes successfully
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
- `results/test/split/` contains:
  `test_rna_Normal_R1.fq.gz`
  `test_rna_Normal_R2.fq.gz`
  `test_rna_Co2_R1.fq.gz`
  `test_rna_Co2_R2.fq.gz`
  `SAM_RG_Header_test_rna_Normal.tsv`
  `SAM_RG_Header_test_rna_Co2.tsv`
- `results/test/usam/` contains:
  `test_rna_Normal_tagged.usam`
  `test_rna_Co2_tagged.usam`
- `results/test_align_mock/align/` contains:
  `test_rna_Normal.Solo.outGeneFull/`
  `test_rna_Normal.filtered_cells.bam`
  `test_rna_Co2.Solo.outGeneFull/`
  `test_rna_Co2.filtered_cells.bam`
- `results/test/pipeline_info/` contains the configured Nextflow report artifacts

The external/local real-RNA validation path is accepted when:

- `nextflow run . -profile test_real_rna --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml --outdir results/test_real_rna` completes successfully
- the real-data files remain external/local and are not committed to this repo
- `results/test_real_rna/tagging/` contains SB, UM, CB, and trim outputs for the supplied sample
- `results/test_real_rna/split/` contains one FASTQ pair plus one `SAM_RG_Header_*.tsv` per RNA SB group in the supplied `rna_sb_group_map`
- `results/test_real_rna/usam/` contains one `<sample>_<group>_tagged.usam` per RNA SB group in the supplied `rna_sb_group_map`
- `results/test_real_rna/pipeline_info/` contains the configured Nextflow report artifacts
- `nextflow run . -profile test_real_rna --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml --outdir results/test_real_rna_align --run_align_rna true --rna_ref_base_dir /mnt/dataFast/ahrmad --rna_align_species human` completes successfully
- `results/test_real_rna_align/align/` contains one `<sample>_<group>.Solo.outGeneFull/` plus one `<sample>_<group>.filtered_cells.bam` per grouped `.usam`

Exact test command:

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

## Troubleshooting

To verify the host prerequisites enforced for every pipeline run:

```bash
bin/check_codon_seq_host.sh
```

That check verifies:

- `codon` is on `PATH`
- Codon reports exactly `0.16.3`
- Seq plugin metadata exists where Codon expects it by default
- Seq plugin metadata reports exactly `0.11.3`
