# AGENTS

- Use Nextflow DSL2 and an nf-core-style layout (`modules/`, `subworkflows/`, `conf/`, `assets/`, `bin/`, `workflows/` if appropriate).
- Keep pipeline runnable end-to-end with a minimal test dataset and `-profile test`.
- Pipeline input is a single YAML samplesheet file (`params.samplesheet`).
- Put HTML reports (`report`, `timeline`, `trace`, `DAG`) under `${params.outdir}/pipeline_info/`.
- Prefer small, reviewable PRs and include acceptance criteria in PR description.
- Add/keep README instructions aligned with actual CLI params.

## Repo-specific requirements

- Preserve current script behavior unless a change is required for pipeline integration.
- Prefer thin Nextflow wrappers around existing scripts in `modules/` and keep script business logic in place.
- Document each module's expected inputs, outputs, and command invocation in comments or module docs.
- Add a minimal YAML samplesheet example under `assets/`.
- Add a minimal test dataset under `assets/testdata/` or `tests/data/`.
- The pipeline is done when `nextflow run . -profile test --samplesheet <test-yaml> --outdir <outdir>` completes successfully and produces expected outputs.
- Keep README command examples synchronized with actual params in `nextflow.config` and workflow code.
- When uncertain about script ordering or file contracts, inspect the code and state assumptions in PR notes.
