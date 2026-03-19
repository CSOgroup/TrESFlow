/*
 * Module: SPLIT_DUPLICATES_DNA
 * Upstream reference:
 *   samtools view --threads <threads> --bam --with-header \
 *     --require-flags 0x400 \
 *     --output <temp_dup.bam> \
 *     --unoutput <sample>_NoDup.bam \
 *     <sample>_MarkedDup.bam
 *
 *   samtools index --threads <threads> --bai \
 *     --output <sample>_NoDup.bam.bai \
 *     <sample>_NoDup.bam
 *
 * Inputs:
 *   - duplicate-marked DNA BAM from GATK MarkDuplicates
 * Outputs:
 *   - duplicate-filtered NoDup BAM
 *   - NoDup BAM index
 *
 * Notes:
 *   - This wrapper preserves the launcher's post-MarkDuplicates duplicate split and NoDup indexing.
 *   - The launcher's temporary DUP.bam sideproduct is created inside the task work dir and removed before publish.
 */

process SPLIT_DUPLICATES_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/dna_nodup", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(markedDupBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam.bai"), emit: bai

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock nodup bam for %s\n' "${splitName}" > "${splitName}_NoDup.bam"
        printf 'mock nodup bai for %s\n' "${splitName}" > "${splitName}_NoDup.bam.bai"
        """
    }
    else {
        """
        if [[ ! -x "${params.runtime_samtools}" ]]; then
          echo "Missing configured DNA runtime executable: ${params.runtime_samtools}" >&2
          exit 1
        fi

        echo "Using SAMTOOLS_BIN=${params.runtime_samtools}"

        "${params.runtime_samtools}" view \\
          --threads "${task.cpus}" \\
          --bam \\
          --with-header \\
          --require-flags 0x400 \\
          --output "${splitName}.DUP.bam" \\
          --unoutput "${splitName}_NoDup.bam" \\
          "${markedDupBam}"

        "${params.runtime_samtools}" index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}_NoDup.bam.bai" \\
          "${splitName}_NoDup.bam"

        rm -f "${splitName}.DUP.bam"
        """
    }
}
