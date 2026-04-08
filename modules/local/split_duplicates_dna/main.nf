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

    input:
    tuple val(splitName), val(meta), path(markedDupBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam.bai"), emit: bai
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock nodup bam for %s\n' "${splitName}" > "${splitName}_NoDup.bam"
        printf 'mock nodup bai for %s\n' "${splitName}" > "${splitName}_NoDup.bam.bai"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        if [[ ! -x "\$SAMTOOLS_BIN" ]]; then
          echo "Missing configured DNA runtime executable: \$SAMTOOLS_BIN" >&2
          exit 1
        fi

        echo "Using SAMTOOLS_BIN=\$SAMTOOLS_BIN"

        "\$SAMTOOLS_BIN" view \\
          --threads "${task.cpus}" \\
          --bam \\
          --with-header \\
          --require-flags 0x400 \\
          --output "${splitName}.DUP.bam" \\
          --unoutput "${splitName}_NoDup.bam" \\
          "${markedDupBam}"

        "\$SAMTOOLS_BIN" index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}_NoDup.bam.bai" \\
          "${splitName}_NoDup.bam"

        rm -f "${splitName}.DUP.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
