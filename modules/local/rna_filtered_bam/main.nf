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

process RNA_FILTERED_BAM {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/align", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(soloDir), path(alignedBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}.filtered_cells.bam"), emit: filtered_bam

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock filtered bam for %s\n' "${splitName}" > "${splitName}.filtered_cells.bam"
        """
    }
    else {
        """
        if [[ ! -x "\$SAMTOOLS_BIN" ]]; then
          echo "Missing configured RNA runtime executable: \$SAMTOOLS_BIN" >&2
          exit 1
        fi

        bash "${params.core_scripts_dir}/RNA_FILTERED_BAM.sh" \\
          "${splitName}" \\
          "${soloDir}" \\
          "${alignedBam}" \\
          "." \\
          "${task.cpus}"
        """
    }
}
