# Phase 0 production regression suite

This suite freezes biological and published-output behavior from the exact
`v1.1.1` tag (`40a383e80d945618952f8e5bfddb73ed3dc63af6`). It is separate from the
fast mock-driven nf-tests: every process in these scenarios executes its real
command.

The fixture generator creates an artificial CC0 genome, annotation, barcode
reads, genomic reads, and sub-megabyte STAR/BWA-MEM2 indices in a caller-owned
temporary directory. No production reference or source data is copied into the
repository. `asset-manifest.json` records every generated file's size and
SHA-256 plus the generator hash, exact index-builder versions, and complete
commands. Fixture generation fails on an index-builder version mismatch.

Three scenarios are covered: RNA only, DNA single tagmentation, and DNA dual
tagmentation. The dual fixture contains one exact linker-artifact pair and all
DNA reads share genomic coordinates so artifact rejection and duplicate
handling are nontrivial.

Run and compare all scenarios with an empty workspace:

```bash
python tests/regression/run_regression.py \
  --workspace /tmp/tresflow-phase0-regression \
  --env-prefix /path/to/v1.1.1/environment \
  --nextflow /path/to/nextflow \
  --python /path/to/python \
  --star /path/to/STAR \
  --bwa-mem2 /path/to/bwa-mem2 \
  --samtools /path/to/samtools
```

Use `--engine docker` or `--engine apptainer` with a separate empty workspace
to exercise declared task containers. At the Phase 0 boundary most custom
processes still have no container directive, so these modes are intentionally
available for migration verification but do not yet prove full isolation from
the host environment.

Golden capture is intentionally separate and overwrite-protected. The initial
capture must add `--capture-baseline`; capture fails if a scenario JSON already
exists. Review its provenance and semantic diff before adding a new baseline.

Normalization rules:

- compare complete relative file and directory structure;
- parse metric TSV/CSV and STAR Matrix Market/barcode/feature/summary files;
- decode BAM with samtools, normalize volatile `@PG` fields and STAR's
  task-allocated `--runThreadN`, sort records and tags, and retain alignment
  fields, sequences, qualities, flags, coordinates, CIGARs, and tags;
- decode bigWig chromosomes, intervals, and values with pyBigWig;
- extract visible TrES report data while excluding CSS and JavaScript;
- compare normalized MultiQC general statistics, raw module data, and exported
  tables rather than its generated HTML/parquet bytes;
- extract FastQC data and summaries from ZIP archives rather than ZIP bytes;
- decompress and sort tag-record gzip content, ignoring gzip headers;
- replace timestamps and absolute paths, and ignore content of Nextflow's
  runtime report/timeline/trace/flowchart while still requiring their paths.

Byte comparison is not used for BAM, bigWig, HTML, MultiQC, ZIP, gzip, or other
tool-generated binary formats.
