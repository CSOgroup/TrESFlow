/*
 * Convert nf-core/samtools/quickcheck exit codes into a small report artifact.
 */

process SAMTOOLS_QUICKCHECK_REPORT {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir ?: "${projectDir}/results"}/qc/samtools", mode: params.publish_dir_mode, overwrite: true, pattern: "*.quickcheck.tsv"

    input:
    tuple val(meta), path(bam), val(exitCode)

    output:
    tuple val(meta), path("*.quickcheck.tsv"), emit: report
    path("versions.yml"), emit: versions

    script:
    def prefix = meta.id
    """
    printf 'id\tbam\texit_code\tstatus\n' > "${prefix}.quickcheck.tsv"
    if [[ "${exitCode}" == "0" ]]; then
      status="pass"
    else
      status="fail"
    fi
    printf '%s\t%s\t%s\t%s\n' "${prefix}" "${bam}" "${exitCode}" "\${status}" >> "${prefix}.quickcheck.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
