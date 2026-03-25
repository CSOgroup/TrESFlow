# AGENTS

- Use Nextflow DSL2 and keep the repo organized with a clear pipeline layout (`modules/`, `subworkflows/`, `conf/`, `assets/`, `bin/`, `workflows/` where appropriate).
- Pipeline input is a single YAML samplesheet file (`params.samplesheet`).
- Put HTML reports (`report`, `timeline`, `trace`, `DAG`) under `${params.outdir}/pipeline_info/`.
- Keep README instructions aligned with the actual supported params and workflow behavior.
- Prefer changes that make the pipeline easier to understand, explain, and maintain.

## Repo-specific requirements

- Treat this repo as the source of truth for the implemented pipeline.
- Preserve validated scientific behavior and output contracts unless an intentional contract change is being made.
- Prefer repo-owned implementations over thin wrappers when that improves clarity and maintainability.
- Backward compatibility is optional unless explicitly requested; lean structure is preferred over legacy support.
- Keep one clear supported public contract for inputs and runtime; avoid parallel legacy interfaces.
- Keep one fast smoke-test path if practical, but do not preserve obsolete mock or compatibility layers just to keep old structure alive.
- When refactoring, optimize for:
  - clarity
  - coherence
  - repo ownership
  - minimal redundancy
- It is acceptable to:
  - remove stale files
  - remove unused configs/profiles
  - simplify docs aggressively
  - collapse compatibility shims
  - internalize code previously inherited from upstream
- Document each module’s expected inputs, outputs, and command invocation in comments or module docs when useful, but do not preserve boilerplate documentation that no longer helps readers.
- When uncertain about workflow ordering or file contracts, inspect the code and state assumptions clearly.
- The pipeline is considered healthy when the currently supported smoke and real-data validation paths complete successfully and produce the expected outputs for the supported contract.
- When there is a tradeoff between backward compatibility and a cleaner supported contract, prefer the cleaner supported contract unless told otherwise.
