/*
 * Delete only explicitly enumerated producer-owned FASTQs after lifecycle
 * completion signals from every enabled downstream consumer have arrived.
 * Paths are values, not staged inputs, so this task cannot create another
 * FASTQ copy while reclaiming the producer's file.
 */

include { runtimeShellExports; intermediateFastqCleanupCommand } from '../runtime_support/main'

process CLEANUP_INTERMEDIATE_FASTQS {
    cache false

    tag "${cleanupId}"
    label 'process_single'

    input:
    tuple val(cleanupId), val(meta), val(fastqs), val(workDir)

    output:
    tuple val(cleanupId), path('fastq_cleanup.tsv'), emit: report

    script:
    def runtimeExports = runtimeShellExports(meta)
    def cleanupCommand = intermediateFastqCleanupCommand(
        true,
        'producer-output',
        workDir,
        "${projectDir}/bin/cleanup_intermediate_fastqs.py",
        fastqs,
        'fastq_cleanup.tsv'
    )

    """
    ${runtimeExports}

    ${cleanupCommand}
    """
}
