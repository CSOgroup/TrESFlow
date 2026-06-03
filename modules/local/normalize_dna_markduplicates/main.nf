/*
 * Normalize nf-core/gatk4/markduplicates outputs back to the TrESFlow DNA
 * filename contract used by downstream NoDup extraction and reporting.
 */

process NORMALIZE_DNA_MARKDUPLICATES {
    tag "${splitName}"
    label 'process_single'

    publishDir "${params.outdir ?: "${projectDir}/results"}/dna_align", mode: params.publish_dir_mode, overwrite: true, pattern: "${splitName}_MarkedDup.bam*"
    publishDir "${params.outdir ?: "${projectDir}/results"}/dna_align", mode: params.publish_dir_mode, overwrite: true, pattern: "${splitName}.DuplicateMetrics.txt"

    input:
    tuple val(splitName), val(meta), path(markedDupBam), path(markedDupBai), path(markedDupMetrics)

    output:
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}.DuplicateMetrics.txt"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    """
    copy_if_needed() {
      src="\$1"
      dest="\$2"
      if [[ "\$(readlink -f "\${src}")" != "\$(readlink -f "\${dest}" 2>/dev/null || true)" ]]; then
        cp -L "\${src}" "\${dest}"
      fi
    }

    copy_if_needed "${markedDupBam}" "${splitName}_MarkedDup.bam"
    copy_if_needed "${markedDupBai}" "${splitName}_MarkedDup.bam.bai"
    copy_if_needed "${markedDupMetrics}" "${splitName}.DuplicateMetrics.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
