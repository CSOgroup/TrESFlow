/*
 * Module: RUN_SC_PROCESS
 * Purpose:
 *   - Run the single optional downstream sc_process.py analysis after the validated
 *     RNA and DNA branches have both finished and the shared staged workdir exists.
 *   - Keep the core pipeline unchanged by requiring the explicit --run_sc_process toggle.
 *
 * Inputs:
 *   - shared stage label
 *   - shared species/genome labels
 *   - staged flat shared workdir produced by STAGE_SC_PROCESS_INPUTS
 *
 * Outputs:
 *   - sc_process_run/ containing the staged inputs plus any downstream outputs emitted
 *     by one single upstream sc_process.py invocation
 *
 * Invocation:
 *   - Copies the staged shared workdir locally
 *   - Exports SNAP_DATA_DIR explicitly to the configured local SnapATAC cache
 *   - Uses task-local HOME, MPLCONFIGDIR, and NUMBA_CACHE_DIR so the optional
 *     downstream runtime stays deterministic and writable inside Nextflow work dirs
 *   - Runs the configured optional sc_process.py entrypoint exactly once
 */

process RUN_SC_PROCESS {
    tag "${stageLabel}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/shared_stage", mode: 'copy', overwrite: true

    input:
    tuple val(stageLabel), val(species), val(genome), path(stageDir)

    output:
    path("sc_process_run"), emit: run_dir

    script:
    """
    if [[ ! -x "\$PYTHON3_BIN" ]]; then
      echo "Missing configured shared runtime executable: \$PYTHON3_BIN" >&2
      exit 1
    fi

    if [[ ! -d "${params.runtime_snap_data_dir}" ]]; then
      echo "Missing configured SnapATAC cache directory: ${params.runtime_snap_data_dir}" >&2
      exit 1
    fi

    if [[ ! -f "${params.optional_sc_process_script}" ]]; then
      echo "Missing configured optional sc_process.py entrypoint: ${params.optional_sc_process_script}" >&2
      exit 1
    fi

    mkdir -p "sc_process_run" \\
             "sc_process_runtime/home" \\
             "sc_process_runtime/mplconfig" \\
             "sc_process_runtime/numba"
    cp -a "${stageDir}/." "sc_process_run/"

    export SNAP_DATA_DIR="${params.runtime_snap_data_dir}"
    export HOME="\$PWD/sc_process_runtime/home"
    export MPLCONFIGDIR="\$PWD/sc_process_runtime/mplconfig"
    export NUMBA_CACHE_DIR="\$PWD/sc_process_runtime/numba"

    {
      printf 'key\\tvalue\\n'
      printf 'SNAP_DATA_DIR\\t%s\\n' "\$SNAP_DATA_DIR"
      printf 'HOME\\t%s\\n' "\$HOME"
      printf 'MPLCONFIGDIR\\t%s\\n' "\$MPLCONFIGDIR"
      printf 'NUMBA_CACHE_DIR\\t%s\\n' "\$NUMBA_CACHE_DIR"
    } > "sc_process_run/runtime_env.tsv"

    echo "Using SNAP_DATA_DIR=\$SNAP_DATA_DIR"
    echo "Using HOME=\$HOME"
    echo "Using MPLCONFIGDIR=\$MPLCONFIGDIR"
    echo "Using NUMBA_CACHE_DIR=\$NUMBA_CACHE_DIR"

    "\$PYTHON3_BIN" "${params.optional_sc_process_script}" \\
      --workdir "sc_process_run" \\
      --mo-map "sc_process_run/mo_map.tsv" \\
      --sb-map-dna "sc_process_run/sb_group_map.tsv" \\
      --sb-map-rna "sc_process_run/sb_group_map.tsv" \\
      --pairs-tsv "sc_process_run/pairs.tsv" \\
      --genome "${genome}" \\
      --species "${species}" \\
      --threads "${task.cpus}"
    """
}
