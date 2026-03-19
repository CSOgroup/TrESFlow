# TrESFlow

TrESFlow is a Nextflow DSL2 wrapper around the read-only source material in `upstream/source_scripts/`. The current implementation is intentionally small: it parses a single YAML samplesheet, runs the implemented RNA workflow through `AlignRNA.sh`, and now also runs the implemented DNA workflow through DNA `MarkDuplicates`.

## Current slice

- `TAG_RNA_SAMPLE_BARCODE` wraps `upstream/source_scripts/Tag.codon` for the RNA sample-barcode step.
- `TAG_RNA_UMI` wraps `upstream/source_scripts/Tag_UMI.codon` for the RNA UMI step.
- `TAG_RNA_CELL_BARCODE` wraps `upstream/source_scripts/Tag_Lig3.codon` for the RNA cell-barcode step.
- `TRIM_RNA_FASTQS` wraps the immediate upstream RNA `trim_galore` step after CB tagging.
- `SPLIT_RNA_READS` wraps the immediate upstream RNA `Split_ReadsV2.codon` step in `rna` mode after trimming.
- `FQ_TO_SAM` wraps the immediate upstream RNA `FqToSAM.codon` step after read splitting.
- `ALIGN_RNA` wraps the immediate upstream RNA `AlignRNA.sh` step after grouped `.usam` generation.
- `TAG_DNA_SAMPLE_BARCODE` wraps the upstream DNA sample-barcode `Tag.codon` step on `I2`.
- `TAG_DNA_MODALITY_BARCODE` wraps the upstream DNA modality-barcode `Tag.codon` step on `I2`.
- `TAG_DNA_CELL_BARCODE` wraps the upstream DNA `Tag_Lig3.codon` ligation-barcode step on `I1`.
- `TRIM_DNA_FASTQS` wraps the immediate upstream DNA `trim_galore` step after CB tagging.
- `SPLIT_DNA_READS` wraps the immediate upstream DNA `Split_ReadsV2.codon` step in `dna` mode after trimming.
- `ALIGN_DNA` wraps the immediate upstream DNA `AlignDNA.sh` step after grouped DNA split outputs.
- `MARK_DUPLICATES_DNA` wraps the immediate upstream DNA `gatk MarkDuplicates` step after direct `AlignDNA.sh` BAM outputs.
- `-profile test` uses lightweight mock wrappers for the wrapped RNA steps, but the pipeline still enforces host Codon `0.16.3` and Seq `0.11.3` before any run starts.
- `-profile test_dna` uses lightweight mock wrappers for the wrapped DNA steps through duplicate-marked DNA BAM outputs, and still enforces host Codon `0.16.3` and Seq `0.11.3` before any run starts.
- Every pipeline run requires host-installed Codon `0.16.3` and Seq `0.11.3`.
- For normal execution on this server, the pipeline binds `python3`, `trim_galore`, `STAR`, `samtools`, `bedGraphToBigWig`, and `bwa-mem2` explicitly to `/home/annan/micromamba/envs/tres/bin/...`.
- `gatk` remains explicitly pinned outside that env at `/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk`.
- `PATH` still prepends `/home/annan/micromamba/envs/tres/bin` as a compatibility fallback for unbound tools and future steps. Override with `--runtime_env_prefix` if needed.
- `envs/first_slice.yml` is the source of truth for software requirements around the current implemented wrappers and mock path.
- A separate `test_real_rna` profile is available for external/local validation data through the full implemented RNA boundary.
- An optional `docker` profile containerizes the current Python wrapper steps for the smoke-test path only.
- Pinning the real-mode toolchain does not change the mock `-profile test` or `-profile test,docker` paths.
- The first true shared downstream step is still one future `sc_process.py` call. This repo does not reproduce the second erroneous `sc_process.py` invocation from `MAINLAUNCH.sh`.
- Downstream RNA `.ubam` conversion, DNA duplicate splitting to `_NoDup.bam`, DNA coverage generation, and shared `sc_process.py` analysis remain intentionally out of scope for the current implementation.
- Every run also writes `pipeline_info/runtime_contract.tsv`, which records the configured explicit runtime binaries plus the host Codon/Seq preflight output.

## Implemented Boundaries

### RNA

The implemented RNA workflow is:

1. `Tag.codon`
2. `Tag_UMI.codon`
3. `Tag_Lig3.codon`
4. `trim_galore`
5. `Split_ReadsV2.codon` in `rna` mode
6. `FqToSAM.codon`
7. `AlignRNA.sh`

This repo implements RNA steps 1 through 7 as the default RNA path.

### DNA

The upstream DNA order relevant to this repo is:

1. `Tag.codon` on `I2` for sample barcode (`SB`)
2. `Tag.codon` on `I2` for modality barcode (`MO`)
3. `Tag_Lig3.codon` on `I1` for ligation/cell barcode (`CB`) and `RG`
4. `trim_galore`
5. `Split_ReadsV2.codon` in `dna` mode
6. `AlignDNA.sh`
7. `gatk MarkDuplicates`
8. downstream duplicate splitting and coverage generation

This repo currently implements the next DNA boundary: steps 1 through 7.
DNA duplicate splitting to `_NoDup.bam`, coverage generation, and shared downstream analysis remain intentionally out of scope.

## Layout

- `main.nf`
- `nextflow.config`
- `conf/base.config`
- `conf/test.config`
- `conf/test_dna.config`
- `workflows/treseq.nf`
- `subworkflows/local/initial_rna_tagging.nf`
- `subworkflows/local/initial_dna_tagging.nf`
- `modules/local/tag_rna_sb/main.nf`
- `modules/local/tag_rna_umi/main.nf`
- `modules/local/tag_rna_cell_barcode/main.nf`
- `modules/local/trim_rna_fastqs/main.nf`
- `modules/local/split_rna_reads/main.nf`
- `modules/local/fq_to_sam/main.nf`
- `modules/local/align_rna/main.nf`
- `modules/local/tag_dna_sb/main.nf`
- `modules/local/tag_dna_modality/main.nf`
- `modules/local/tag_dna_cell_barcode/main.nf`
- `modules/local/trim_dna_fastqs/main.nf`
- `modules/local/split_dna_reads/main.nf`
- `modules/local/align_dna/main.nf`
- `modules/local/mark_duplicates_dna/main.nf`
- `assets/samplesheet.example.yaml`
- `assets/samplesheet.dna.example.yaml`
- `assets/samplesheet.real_rna.template.yaml`
- `assets/testdata/`
- `scripts/install_codon_0.16.3.sh`

## Required CLI params

- `--samplesheet`
- `--outdir`

If the samplesheet contains any RNA samples, real and mock RNA runs also require:

- `--rna_ref_base_dir`
  Must contain either `GRCh38_TrES/star` plus `hg38.chrom.sizes`, or `GRCm39_TrES/star` plus `mm39.chrom.sizes`, depending on species.
- `--rna_align_species`
  Supported values: `human`, `mouse`

If the samplesheet contains any DNA samples, real and mock DNA runs also require:

- `--dna_bwa_reference`
  bwa-mem2 index prefix. The path itself may be a prefix rather than a regular file; the required sidecars are `${prefix}.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, and `.pac`.
- `--dna_blacklist_bed`
  BED file passed directly to `AlignDNA.sh` for blacklist filtering.
- `--dna_effective_genome_size`
  Integer passed directly to `AlignDNA.sh`.

## Optional CLI params

- `--upstream_dir`
  Default: `./upstream/source_scripts`
- `--runtime_env_prefix`
  Default: `/home/annan/micromamba/envs/tres`
  If `${runtime_env_prefix}/bin` exists, every task prepends it to `PATH` before running.
- `--runtime_python`
  Default: `/home/annan/micromamba/envs/tres/bin/python3`
- `--runtime_trim_galore`
  Default: `/home/annan/micromamba/envs/tres/bin/trim_galore`
- `--runtime_star`
  Default: `/home/annan/micromamba/envs/tres/bin/STAR`
- `--runtime_samtools`
  Default: `/home/annan/micromamba/envs/tres/bin/samtools`
- `--runtime_bedgraph_to_bigwig`
  Default: `/home/annan/micromamba/envs/tres/bin/bedGraphToBigWig`
- `--runtime_bwa_mem2`
  Default: `/home/annan/micromamba/envs/tres/bin/bwa-mem2`
- `--runtime_bam_coverage`
  Default: `/home/annan/micromamba/envs/tres/bin/bamCoverage`
  Reserved for the future DNA coverage step; it is recorded in `pipeline_info/runtime_contract.tsv` but not used by the current implemented boundary.
- `--gatk_root`
  Default: `/mnt/dataFast/ahrmad/gatk-4.6.0.0`
  `MARK_DUPLICATES_DNA` runs `${gatk_root}/gatk MarkDuplicates`.
- `--max_cpus`
  Default: `40`
  Total CPU budget for the local executor. Real runs distribute that budget as `ALIGN_RNA=int(max_cpus/2)`, `TRIM_RNA_FASTQS=min(8, int(max_cpus/5))`, `TRIM_DNA_FASTQS=min(8, int(max_cpus/5))`, `SPLIT_RNA_READS=min(4, int(max_cpus/10))`, `SPLIT_DNA_READS=min(4, int(max_cpus/10))`, `ALIGN_DNA=max_cpus`, `MARK_DUPLICATES_DNA=1`, and `1` CPU for the remaining wrapped processes.
  `ALIGN_DNA` takes the full executor budget so the local executor runs one DNA align task at a time; this keeps the wrapper honest because upstream `AlignDNA.sh` hardcodes its own internal thread counts.

## Samplesheet schema

The parser accepts `modality: rna` and `modality: dna` rows in the same top-level YAML file.
Top-level `sb_group_map` is required whenever the samplesheet contains RNA or DNA samples, because both implemented paths use launcher-style sample-barcode grouping for `Split_ReadsV2.codon`.
Top-level `dna_mo_map` is required whenever the samplesheet contains DNA samples.

### RNA example

```yaml
library_name: run_library_name
sb_group_map: path/to/sb_group_map.tsv

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

The shared sample-barcode group map is the launcher-style TSV used by `Split_ReadsV2.codon`: `sample<TAB>sb_group<TAB>sb_bc`.
For the current RNA and DNA workflows it is also the single source of truth for the sample-barcode whitelist passed into `Tag.codon`.
Blank lines, `#` comments, and a literal `sample sb_group sb_bc` header row are ignored.
The pipeline fails if the map has no rows for a sample or if the same `sb_bc` is assigned to multiple groups for that sample.

### DNA example

```yaml
library_name: run_library_name
sb_group_map: path/to/sb_group_map.tsv
dna_mo_map: path/to/mo_map.tsv

samples:
  - id: sample_id
    modality: dna
    reads:
      i1: path/to/I1.fastq.gz
      i2: path/to/I2.fastq.gz
      r1: path/to/R1.fastq.gz
      r2: path/to/R2.fastq.gz
    barcodes:
      sample:
        bc_len: 4
        bc_start: 14
        hd: 1
        tag: SB
        first_pass: first_pass
        reverse_complement: true
      modality:
        whitelist: path/to/dna_modality_whitelist.txt
        bc_len: 8
        bc_start: 18
        hd: 1
        tag: MO
        first_pass: not_first_pass
        reverse_complement: true
      cell:
        whitelist: path/to/ligation_whitelist.txt
        bc_len: 8
        hd: 1
        tag: CB
```

The current DNA slice follows the launcher-supported barcode sources and order exactly:
sample and modality barcodes come from `I2`, then ligation/cell barcodes come from `I1`.
For DNA sample-barcode tagging and DNA split, `sb_group_map` is the shared sample-barcode grouping TSV and the single source of truth for experiment-used DNA sample barcodes.
`dna_mo_map` is the launcher-style modality map TSV consumed by `Split_ReadsV2.codon` in `dna` mode.
The supported DNA MO-map form for this repo is the launcher-style 4-column TSV: `sample<TAB>sb_group<TAB>mark<TAB>mo_bc`.
DNA alignment uses explicit CLI inputs matching `AlignDNA.sh`: `--dna_bwa_reference`, `--dna_blacklist_bed`, and `--dna_effective_genome_size`.
DNA duplicate-marking, coverage, and downstream QC remain intentionally out of scope.

## Running the test profile

```bash
nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test
```

Expected tagging and trimming outputs land under `results/test/tagging/`, split outputs land under `results/test/split/`, grouped unmapped SAM outputs land under `results/test/usam/`, mocked alignment outputs land under `results/test/align/`, and Nextflow report artifacts land under `results/test/pipeline_info/`. The trace artifact is the standard Nextflow tabular trace file, written as `execution_trace.tsv` in that directory.
The upstream launcher deletes the untrimmed CB FASTQs after `trim_galore`; this slice keeps them published as intermediates and advances on the trimmed `_val_1` / `_val_2` outputs.
`Split_ReadsV2.codon` has one ambiguity in its comments versus examples: the code comments discuss dropping an injected leading base from `SB`, but the upstream RNA map example uses full `SB` strings. This repo follows the script's actual lookup behavior: raw `SB` match first, then drop-first fallback.

## Running the DNA test profile

```bash
nextflow run . -profile test_dna --samplesheet assets/samplesheet.dna.example.yaml --outdir results/test_dna
```

Expected DNA outputs land under `results/test_dna/dna_tagging/`, `results/test_dna/dna_split/`, and `results/test_dna/dna_align/`:
sample-barcode-tagged FASTQs, sample-plus-modality-tagged FASTQs, SB/MO/CB count and stats files, CB tag records, trimmed `_val_1` / `_val_2` DNA FASTQs, split per-group per-mark FASTQ pairs, `SAM_RG_Header_*.tsv` files, mocked aligned BAMs, mocked BAM indexes, and mocked per-barcode count TSVs.

## Real RNA Validation

The real-data validation path is intentionally external/local. No real RNA fixture is committed to this repo.

Minimum real-input contract for the implemented RNA workflow:

- top-level `library_name`
- top-level `sb_group_map`
- per-sample `id` matching the `sample` column in `sb_group_map` for the groups you want emitted.
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

Use `assets/samplesheet.real_rna.template.yaml` as the committed template. The real files above are expected to remain external/local.
The template keeps `samples[0].id: day15` because the sample id must match the first `sample` column in `sb_group_map`.
The current external/local example still uses the legacy on-disk filename `sb_map_RNA.tsv`, but the pipeline contract treats it as a shared sample-barcode group map.

Exact run command for host-based real-RNA validation:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna \
  --rna_ref_base_dir /path/to/reference_base \
  --rna_align_species human
```

Optional variant if you want Nextflow to provision `trim_galore` and Python task dependencies through the checked-in conda environment while still using the globally required host Codon `0.16.3` plus Seq `0.11.3`:

```bash
nextflow run . \
  -profile test_real_rna,conda_dev \
  --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna \
  --rna_ref_base_dir /path/to/reference_base \
  --rna_align_species human
```

Expected outputs through the implemented RNA boundary:

- `outdir/tagging/` with SB, UM, CB, and trim outputs
- `outdir/split/` with one FASTQ pair and one `SAM_RG_Header_*.tsv` per RNA SB group in `sb_group_map`
- `outdir/usam/` with one `<sample>_<group>_tagged.usam` per RNA SB group in `sb_group_map`
- `outdir/align/` with one `<sample>_<group>.Solo.outGeneFull/` STARsolo directory per grouped `.usam`
- `outdir/align/` with one `<sample>_<group>.filtered_cells.bam` per grouped `.usam`
- `outdir/align/` with stranded and unstranded bigWig tracks per grouped `.usam` when `AlignRNA.sh` emits them
- `outdir/pipeline_info/` with the configured Nextflow report artifacts
- `outdir/pipeline_info/runtime_contract.tsv` with the configured explicit runtime binaries and the host Codon/Seq preflight details

## Real DNA Validation

The real DNA validation path for this pass uses the provided real DNA inputs under `assets/test_realdata/` and stops at the direct `gatk MarkDuplicates` output boundary.
The pipeline requires the shared sample-barcode group map plus the DNA modality map for DNA split, and explicit alignment CLI inputs matching `AlignDNA.sh`.
For DNA sample-barcode tagging, `sb_group_map` is also the single source of truth for the effective SB whitelist passed into upstream `Tag.codon`.

Exact run command for host-based real-DNA validation:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.dna.RealDATAexample.yaml \
  --outdir results/test_dna_real_markdup \
  --dna_bwa_reference /path/to/bwa_index_prefix \
  --dna_blacklist_bed /path/to/blacklist.bed \
  --dna_effective_genome_size 2913022398
```

The current provided `assets/test_realdata/sb_map_RNA.tsv` filename is legacy, but it is consumed as the generic `sb_group_map`.
On this server, deriving the DNA SB whitelist from `sb_group_map` plus the explicit launcher-style human alignment inputs makes the provided real DNA run complete through `MarkDuplicates`.
`AlignDNA.sh` directly emits the filtered BAM, BAM index, and properly paired mapped reads per barcode TSV; it also hardcodes `min_good_reads_in_cells=100`.
`MARK_DUPLICATES_DNA` preserves the launcher `gatk MarkDuplicates` contract with `--BARCODE_TAG CB`, `--REMOVE_DUPLICATES false`, `--CREATE_INDEX true`, and `--MAX_RECORDS_IN_RAM 10000000`.

Expected outputs through the implemented DNA boundary when `sb_group_map`, `dna_mo_map`, and the explicit DNA alignment inputs satisfy the upstream contract:

- `outdir/dna_tagging/` with one pair of SB-tagged FASTQs
- `outdir/dna_tagging/` with one pair of SB-plus-MO-tagged FASTQs
- `outdir/dna_tagging/` with one pair of SB-plus-MO-plus-CB-tagged FASTQs
- `outdir/dna_tagging/` with `Reads_Per_Barcode`-style counts plus barcode statistics and tag-record TSVs from the three wrapped upstream DNA steps
- `outdir/dna_tagging/` with trimmed `_val_1` / `_val_2` FASTQs
- `outdir/dna_split/` with one `<sample>_<group>_<mark>_R1.fq.gz` and `<sample>_<group>_<mark>_R2.fq.gz` per valid group and mark combination
- `outdir/dna_split/` with one `SAM_RG_Header_<sample>_<group>_<mark>.tsv` per valid group and mark combination
- `outdir/dna_align/` with one `<sample>_<group>_<mark>.bam` per valid group and mark combination
- `outdir/dna_align/` with one `<sample>_<group>_<mark>.bam.bai` per valid group and mark combination
- `outdir/dna_align/` with one `<sample>_<group>_<mark>_ProperPairedMapped_reads_per_barcode.tsv` per valid group and mark combination
- `outdir/dna_dedup/` with one `<sample>_<group>_<mark>_MarkedDup.bam` per valid group and mark combination
- `outdir/dna_dedup/` with one `<sample>_<group>_<mark>_MarkedDup.bam.bai` per valid group and mark combination
- `outdir/dna_dedup/` with one `<sample>_<group>_<mark>.DuplicateMetrics.txt` per valid group and mark combination
- `outdir/pipeline_info/` with the configured Nextflow report artifacts

## Acceptance Criteria

- `nextflow run . -profile test --samplesheet assets/samplesheet.example.yaml --outdir results/test` succeeds through mocked RNA `AlignRNA.sh` outputs.
- `nextflow run . -profile test_dna --samplesheet assets/samplesheet.dna.example.yaml --outdir results/test_dna` succeeds through mocked DNA duplicate-marking outputs.
- `nextflow run . --samplesheet assets/samplesheet.dna.RealDATAexample.yaml --outdir results/test_dna_real_markdup --dna_bwa_reference /path/to/bwa_index_prefix --dna_blacklist_bed /path/to/blacklist.bed --dna_effective_genome_size 2913022398` succeeds through direct `gatk MarkDuplicates` outputs on a host with Codon `0.16.3`, Seq `0.11.3`, `trim_galore`, `bwa-mem2`, `samtools`, and `gatk`.
- The DNA sample-barcode whitelist is derived from `sb_group_map`; a separate DNA sample-barcode whitelist is not part of the pipeline contract.
- RNA output locations remain `tagging/`, `split/`, `usam/`, `align/`, and `pipeline_info/`.
- DNA output locations for the current slice are `dna_tagging/`, `dna_split/`, `dna_align/`, `dna_dedup/`, and `pipeline_info/`.
- The pipeline still writes Nextflow `timeline`, `report`, `trace`, and `DAG` artifacts under `${params.outdir}/pipeline_info/`.
- The pipeline also writes `${params.outdir}/pipeline_info/runtime_contract.tsv` so the configured explicit runtime binaries are easy to verify after a run.

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

### Host + DNA Test

Supported on native Linux and macOS hosts:

```bash
nextflow run . -profile test_dna --samplesheet assets/samplesheet.dna.example.yaml --outdir results/test_dna
```

This path uses the bundled mock wrappers for the current DNA boundary through duplicate-marked DNA outputs, and startup still fails unless host Codon `0.16.3` and Seq `0.11.3` are installed.
This mock path emits `dna_split/`, `dna_align/`, and `dna_dedup/` outputs.

### Host + Real RNA Validation

Supported on native Linux and macOS hosts with external/local RNA inputs:

```bash
nextflow run . \
  -profile test_real_rna \
  --samplesheet /path/to/real_rna_validation/samplesheet.real_rna.yaml \
  --outdir results/test_real_rna \
  --rna_ref_base_dir /path/to/reference_base \
  --rna_align_species human
```

This path runs the implemented RNA workflow in real mode through grouped `AlignRNA.sh` outputs.
It requires host-installed Codon `0.16.3` and Seq `0.11.3`, and by default it executes `python3`, `trim_galore`, `STAR`, `samtools`, and `bedGraphToBigWig` from `/home/annan/micromamba/envs/tres/bin/`.
It expects the real input files to remain external/local and not committed to this repo.

It also requires the reference-base contract expected by `AlignRNA.sh`:
`<ref_base_dir>/GRCh38_TrES/star` with `<ref_base_dir>/hg38.chrom.sizes` for `human`, or `<ref_base_dir>/GRCm39_TrES/star` with `<ref_base_dir>/mm39.chrom.sizes` for `mouse`.

### Host + Real DNA Validation

Supported on native Linux and macOS hosts with the provided real DNA files under `assets/test_realdata/`:

```bash
nextflow run . \
  --samplesheet assets/samplesheet.dna.RealDATAexample.yaml \
  --outdir results/test_dna_real_markdup \
  --dna_bwa_reference /path/to/bwa_index_prefix \
  --dna_blacklist_bed /path/to/blacklist.bed \
  --dna_effective_genome_size 2913022398
```

This path runs the current DNA workflow in real mode through direct `gatk MarkDuplicates` outputs.
It requires host-installed Codon `0.16.3` and Seq `0.11.3`, and by default it executes `python3`, `trim_galore`, `bwa-mem2`, and `samtools` from `/home/annan/micromamba/envs/tres/bin/`.
It uses the top-level `sb_group_map` plus `dna_mo_map` entries from the samplesheet to satisfy the upstream `Split_ReadsV2.codon` DNA contract, and it uses the explicit DNA alignment CLI params to satisfy the upstream `AlignDNA.sh` contract.
For DNA sample-barcode tagging, the pipeline derives the effective SB whitelist directly from `sb_group_map`.
It uses `${params.gatk_root}/gatk` for duplicate marking; on this server the default is `/mnt/dataFast/ahrmad/gatk-4.6.0.0/gatk`.
Duplicate splitting to `_NoDup.bam`, coverage generation, and shared downstream analysis remain intentionally out of scope.

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

This keeps the current wrapper tasks containerized, but it is no longer a standalone portable smoke-test path because startup still enforces host Codon `0.16.3` and Seq `0.11.3`.
Under `-profile docker`, every process labeled `codon_wrapper` runs in the local image `tresflow-first-slice:py312`, built from `docker/first_slice.Dockerfile`.
That currently includes the wrapped RNA path plus the wrapped DNA path through direct `gatk MarkDuplicates` outputs.
This does not make any execution mode fully portable, because Codon `0.16.3` and Seq `0.11.3` are still required on the host in the current implementation.

## Real Mode Host Prerequisites

Every pipeline run currently requires the following on the host:

- native Linux or macOS
- Codon `0.16.3` installed under `${HOME}/.codon` and available on `PATH`
- Seq `0.11.3` installed separately under `${HOME}/.codon/lib/codon/plugins/seq`
- `/home/annan/micromamba/envs/tres/bin/python3`
- `/home/annan/micromamba/envs/tres/bin/trim_galore`
- `/home/annan/micromamba/envs/tres/bin/STAR`
- `/home/annan/micromamba/envs/tres/bin/samtools`
- `/home/annan/micromamba/envs/tres/bin/bedGraphToBigWig`
- `/home/annan/micromamba/envs/tres/bin/bwa-mem2`
- `gatk` available under `${params.gatk_root}/gatk` for DNA duplicate marking runs
- the read-only upstream scripts under `upstream/source_scripts/`

For the `test_real_rna` profile specifically, you also need an external/local samplesheet plus the referenced `I1`, `R1`, `R2`, cell whitelist, and shared sample-barcode group map files. The real-RNA sample id must match the `sample` column in that map. For the current example, that value is `day15`.
For real DNA alignment runs, you also need the referenced `I1`, `I2`, `R1`, `R2`, DNA modality whitelist, ligation whitelist, `sb_group_map`, and `dna_mo_map` files, plus `--dna_bwa_reference`, `--dna_blacklist_bed`, and `--dna_effective_genome_size`.
By default on this server, the pipeline binds the current implemented runtime tools explicitly to `/home/annan/micromamba/envs/tres/bin`, and also prepends that directory to `PATH` as a compatibility fallback. Override the individual `--runtime_*` params if you need different binaries.
`--dna_bwa_reference` is a bwa-mem2 index prefix; the base path itself may be a prefix rather than a regular file, but the sidecars `${prefix}.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, and `.pac` must exist.

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
nextflow run . \
  --samplesheet /path/to/samplesheet.yaml \
  --outdir /path/to/results \
  --rna_ref_base_dir /path/to/reference_base \
  --rna_align_species human
```

For every execution mode, `envs/first_slice.yml` is not enough by itself. It does not install Codon `0.16.3` or Seq `0.11.3`, which remain mandatory host prerequisites until they are containerized.

## Dependency Mapping for the Current Slice

Relevant entries from `envs/first_slice.yml` for the implemented RNA and DNA wrappers:

- Launcher environment when you manually `micromamba activate tres`:
  `nextflow`, `openjdk`
- Current wrapper runtime for the DNA and RNA Python wrapper modules under `modules/local/`:
  explicit `/home/annan/micromamba/envs/tres/bin/python3`
- Current real trimming runtime for `modules/local/trim_rna_fastqs/main.nf` and `modules/local/trim_dna_fastqs/main.nf`:
  explicit `/home/annan/micromamba/envs/tres/bin/trim_galore`, plus its bundled `cutadapt`
- Current real RNA alignment runtime for `modules/local/align_rna/main.nf`:
  explicit `/home/annan/micromamba/envs/tres/bin/STAR`, `/home/annan/micromamba/envs/tres/bin/samtools`, `/home/annan/micromamba/envs/tres/bin/bedGraphToBigWig`

Not provided by `envs/first_slice.yml`:

- Codon `0.16.3`
- Seq `0.11.3` plugin files under `${HOME}/.codon/lib/codon/plugins/seq`

Current process usage:

- `TAG_RNA_SAMPLE_BARCODE` and `TAG_DNA_SAMPLE_BARCODE` run `bin/run_tag.py`.
  Today it only requires Python's standard library in both `mock` and `real` mode.
  In `real` mode it also shells out to host-provided `codon` with `-plugin seq`, which is not currently supplied by `envs/first_slice.yml`.
- `TAG_RNA_UMI` runs [`bin/run_tag_umi.py`](/Users/aannan/GitAA/TrESFlow/bin/run_tag_umi.py).
  It has the same dependency pattern: Python stdlib plus host-provided `codon` with the Seq plugin for `real` mode.
- `TAG_RNA_CELL_BARCODE` and `TAG_DNA_CELL_BARCODE` run `bin/run_tag_lig3.py`.
  It also uses Python stdlib in `mock` mode and shells out to host-provided `codon` with the Seq plugin in `real` mode.
- `TAG_DNA_MODALITY_BARCODE` also runs `bin/run_tag.py` with the same host Codon plus Seq dependency pattern as the sample-barcode wrappers.
- `TRIM_RNA_FASTQS` and `TRIM_DNA_FASTQS` run `bin/run_trim_galore.py`.
  In `mock` mode it gzip-copies the CB-tagged FASTQs to the expected trim_galore `_val_1` / `_val_2` outputs.
  In `real` mode it shells out to the configured explicit `trim_galore` binary with the upstream launcher settings `--quality 10 --gzip --length 20 --paired`.
- `SPLIT_RNA_READS` runs [`bin/run_split_reads_rna.py`](/Users/aannan/GitAA/TrESFlow/bin/run_split_reads_rna.py).
  In `mock` mode it reproduces the upstream RNA grouping behavior from the trimmed FASTQs and writes launcher-style `sample_group_R1/R2.fq.gz` plus `SAM_RG_Header_sample_group.tsv` files.
  In `real` mode it shells out to `codon -plugin seq` for `Split_ReadsV2.codon` in `rna` mode with the launcher-style RNA SB-group map.
- `FQ_TO_SAM` runs [`bin/run_fq_to_sam.py`](/Users/aannan/GitAA/TrESFlow/bin/run_fq_to_sam.py).
  In `mock` mode it reproduces the checked-in `FqToSAM.codon` behavior from split FASTQs and writes grouped unmapped SAM files.
  In `real` mode it shells out to `codon -plugin seq` for `FqToSAM.codon`, which accepts `.fq.gz` inputs directly.
- `ALIGN_RNA` runs [`upstream/source_scripts/AlignRNA.sh`](/Users/aannan/GitAA/TrESFlow/upstream/source_scripts/AlignRNA.sh) directly.
  In `mock` mode it writes placeholder STARsolo directories, filtered BAMs, and bigWig files for smoke-test coverage.
  In `real` mode it binds `STAR`, `samtools`, and `bedGraphToBigWig` explicitly to the configured runtime paths, using the grouped `.usam` plus the reference-base and species contract expected by the upstream shell script.
- The default local CPU budget is `--max_cpus 40`.
  That budget is enforced through the local executor and distributed as:
  `ALIGN_RNA=int(max_cpus/2)`, `TRIM_RNA_FASTQS=min(8, int(max_cpus/5))`, `TRIM_DNA_FASTQS=min(8, int(max_cpus/5))`, `SPLIT_RNA_READS=min(4, int(max_cpus/10))`, and `1` CPU for the remaining wrapped processes.

Current Docker process coverage:

- `TAG_RNA_SAMPLE_BARCODE` is containerized for the smoke-test path
- `TAG_RNA_UMI` is containerized for the smoke-test path
- `TAG_RNA_CELL_BARCODE` is containerized for the smoke-test path
- `TRIM_RNA_FASTQS` is containerized for the smoke-test path
- `SPLIT_RNA_READS` is containerized for the smoke-test path
- `FQ_TO_SAM` is containerized for the smoke-test path
- `ALIGN_RNA` is containerized for the smoke-test path
- `TAG_DNA_SAMPLE_BARCODE` is containerized for the smoke-test path
- `TAG_DNA_MODALITY_BARCODE` is containerized for the smoke-test path
- `TAG_DNA_CELL_BARCODE` is containerized for the smoke-test path
- `TRIM_DNA_FASTQS` is containerized for the smoke-test path
- Codon and Seq are not containerized in this pass

Dependencies present in `envs/first_slice.yml` but currently unused by the implemented boundaries:

- Alignment and downstream genomics tools:
  `samtools`, `star`, `bwa-mem2`, `bedtools`, `fastqc`, `multiqc`, `deeptools`, `ucsc-bedgraphtobigwig`
- Downstream analysis stack for later steps such as `sc_process.py`:
  `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, `matplotlib`, `matplotlib-venn`, `upsetplot`, pip `snapatac2`, pip `MACS3`
- Test and developer tooling not used by the current runtime path:
  `pytest`, `pytest-timeout`
- Other packages currently unused by the implemented boundaries:
  `pyyaml`, `coreutils`, `parallel`, `pigz`, `pbzip2`

## Docker Status

Supported via Docker now:

- An optional `docker` profile exists for the current smoke-test path.
- It containerizes only the Python wrapper processes using the local image `tresflow-first-slice:py312`.
- `docker + test` still containerizes the current wrapper tasks, but it now also requires host Codon `0.16.3` and Seq `0.11.3` because startup enforces the pinned host toolchain globally.

What remains host-dependent right now:

- Codon `0.16.3` for real non-mock execution of the five wrapped upstream Codon RNA steps
- Seq `0.11.3` installation under `${HOME}/.codon/lib/codon/plugins/seq`
- system shell tools such as `bash`, `sort`, `awk`, `grep`, `head`, `tail`, `cat`, `mv`, and `rm`
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
- Real RNA execution is host-only in the current implementation, because the checked-in Docker image does not include the explicitly configured server binaries for `STAR`, `samtools`, or `bedGraphToBigWig`.

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
