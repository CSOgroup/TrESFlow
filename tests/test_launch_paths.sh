#!/usr/bin/env bash

set -euo pipefail

parser="${1:-v2}"
nextflow_bin="${NEXTFLOW_BIN:-nextflow}"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
launch_dir="$(mktemp -d "${TMPDIR:-/tmp}/tresflow-launch-paths.XXXXXX")"

cleanup() {
    if [[ "${TRESFLOW_KEEP_LAUNCH_TMP:-false}" == "true" ]]; then
        printf 'retained launch-path test directory: %s\n' "${launch_dir}" >&2
    else
        rm -rf "${launch_dir}"
    fi
}
trap cleanup EXIT

snapshot_project_results() {
    if [[ -d "${repo_root}/results" ]]; then
        find "${repo_root}/results" -printf '%P\t%y\t%s\t%T@\n' | sort
    fi
}

mkdir -p "${launch_dir}/inputs"
cp "${repo_root}/assets/samplesheet.example.yaml" "${launch_dir}/inputs/samplesheet.yaml"
ln -s "${repo_root}/assets/testdata" "${launch_dir}/inputs/testdata"
cp -a "${repo_root}/scripts/core_runtime" "${launch_dir}/core-runtime"

snapshot_project_results > "${launch_dir}/project-results.before"

config_output="$({
    cd "${launch_dir}"
    env \
        NXF_OFFLINE=true \
        NXF_SYNTAX_PARSER="${parser}" \
        "${nextflow_bin}" config -flat "${repo_root}"
})"
grep -F "timeline.file = '${launch_dir}/results/pipeline_info/execution_timeline.html'" <<< "${config_output}" > /dev/null

(
    cd "${launch_dir}"
    env \
        NXF_OFFLINE=true \
        NXF_SYNTAX_PARSER="${parser}" \
        NXF_WORK="${launch_dir}/work" \
        "${nextflow_bin}" run "${repo_root}/main.nf" \
        -profile test \
        --samplesheet inputs/samplesheet.yaml \
        --outdir relative-results \
        --core_scripts_dir core-runtime \
        --cleanup_work false
)

test -s "${launch_dir}/relative-results/rna_split_fastqs/test_rna_Normal_R1.fastq.gz"
test -s "${launch_dir}/relative-results/pipeline_info/runtime_contract.tsv"
grep -F $'runtime_tmpdir\t'"${launch_dir}/relative-results" \
    "${launch_dir}/relative-results/pipeline_info/runtime_contract.tsv" > /dev/null

# A project-owned wrapper must still be resolved from projectDir, while the
# explicit core runtime override must be resolved from launchDir.
grep -R -F "${repo_root}/bin/run_tag.py" "${launch_dir}/work" --include='.command.sh' > /dev/null
grep -R -F "${launch_dir}/core-runtime/Tag.codon" "${launch_dir}/work" --include='.command.sh' > /dev/null

snapshot_project_results > "${launch_dir}/project-results.after"
cmp "${launch_dir}/project-results.before" "${launch_dir}/project-results.after"

printf 'launch-path regression passed with parser %s\n' "${parser}"
