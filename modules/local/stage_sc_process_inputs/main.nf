/*
 * Module: STAGE_SC_PROCESS_INPUTS
 * Purpose:
 *   - Build the first true shared downstream boundary after the current RNA and DNA branches.
 *   - Stage a flat workdir layout that matches what upstream sc_process.py expects.
 *   - Remain optional so the validated RNA and DNA branches stay useful on their own.
 *
 * Inputs:
 *   - launcher-style DNA modality map TSV
 *   - shared sample-barcode group map TSV
 *   - grouped RNA STARsolo directories
 *   - grouped RNA filtered BAMs
 *   - grouped DNA NoDup BAMs
 *   - shared species/genome labels for readiness notes
 *
 * Outputs:
 *   - flat shared stage directory containing:
 *       * copied RNA Solo.outGeneFull directories
 *       * linked/copied RNA filtered BAMs
 *       * linked/copied DNA NoDup BAMs and BAIs
 *       * mo_map.tsv, sb_group_map.tsv, pairs.tsv, and stage manifests
 *
 * Notes:
 *   - This is intentionally the strongest truthful boundary before one future sc_process.py call.
 *   - The stage dir includes a readiness note describing the current observed blocker
 *     after the SnapATAC cache/runtime bootstrap fix on this server.
 */

import WorkflowSupport

process STAGE_SC_PROCESS_INPUTS {
    tag "${stageLabel}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/shared_stage", mode: 'copy', overwrite: true

    input:
    tuple val(stageLabel), val(species), val(genome), path(moMap), path(sbGroupMap), path(rnaSoloDirs), path(rnaFilteredBams), path(dnaNoDupBams)

    output:
    path("sc_process_stage"), emit: stage_dir

    script:
    def soloArgs = WorkflowSupport.asPathList(rnaSoloDirs).collect { path ->
        "--rna-solo-dir \"${path}\""
    }.join(" \\\n          ")
    def filteredBamArgs = WorkflowSupport.asPathList(rnaFilteredBams).collect { path ->
        "--rna-filtered-bam \"${path}\""
    }.join(" \\\n          ")
    def dnaNoDupArgs = WorkflowSupport.asPathList(dnaNoDupBams).collect { path ->
        "--dna-nodup-bam \"${path}\""
    }.join(" \\\n          ")

    """
    if [[ ! -x "${params.runtime_python}" ]]; then
      echo "Missing configured shared runtime executable: ${params.runtime_python}" >&2
      exit 1
    fi

    "${params.runtime_python}" "${projectDir}/bin/run_stage_sc_process.py" \\
      --stage-dir "sc_process_stage" \\
      --mo-map "${moMap}" \\
      --sb-group-map "${sbGroupMap}" \\
      --species "${species}" \\
      --genome "${genome}" \\
      ${soloArgs} \\
      ${filteredBamArgs} \\
      ${dnaNoDupArgs}
    """
}
