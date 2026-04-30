/*
 * Module: RNA_FILTERED_BAM
 * Runtime command:
 *   bash scripts/core_runtime/RNA_FILTERED_BAM.sh \
 *     <split_name> <solo_dir> <aligned.bam> <outdir> <threads>
 *
 * Inputs:
 *   - STARsolo GeneFull directory from RNA_STARSOLO_ALIGN
 *   - STAR coordinate-sorted aligned BAM from RNA_STARSOLO_ALIGN
 * Outputs:
 *   - filtered-cells RNA BAM
 */

import RuntimeSupport

process RNA_FILTERED_BAM {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/align", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(soloDir), path(alignedBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}.filtered_cells.bam"), emit: filtered_bam
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')
    def runtimeExports = RuntimeSupport.shellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock filtered bam for %s\n' "${splitName}" > "${splitName}.filtered_cells.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        ${runtimeExports}

        if [[ ! -x "\$SAMTOOLS_BIN" ]]; then
          echo "Missing configured RNA runtime executable: \$SAMTOOLS_BIN" >&2
          exit 1
        fi

        bash "${coreScriptsDir}/RNA_FILTERED_BAM.sh" \\
          "${splitName}" \\
          "${soloDir}" \\
          "${alignedBam}" \\
          "." \\
          "${task.cpus}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
