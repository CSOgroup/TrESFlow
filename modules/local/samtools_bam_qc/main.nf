/*
 * Run the TrESFlow BAM sidecar QC suite in one task per BAM. RNA BAMs do not
 * have an index on this channel, so idxstats is emitted only for indexed DNA
 * BAMs. A non-zero quickcheck status is reported without itself failing the
 * task; flagstat, stats, and DNA idxstats retain their normal failure behavior.
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir } from '../runtime_support/main'

process SAMTOOLS_BAM_QC {
    tag "${meta.id}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/TrES_Stats/qc/samtools" }, mode: params.publish_dir_mode, overwrite: true

    input:
    tuple val(meta), path(bam), path(bai), val(runIdxstats)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.idxstats"), optional: true, emit: idxstats
    tuple val(meta), path("*.quickcheck.tsv"), emit: quickcheck
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def runtimeExports = runtimeShellExports(meta)
    def coreScriptsDir = runtimeCoreScriptsDir()

    """
    ${runtimeExports}

    bash "${coreScriptsDir}/SamtoolsBamQc.sh" \\
      "${bam}" \\
      "${prefix}" \\
      "${runIdxstats}" \\
      "${task.cpus}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def idxstatsStub = runIdxstats ? "touch \"${prefix}.idxstats\"" : ''

    """
    printf '%s\\n' \
      '1000000 + 0 in total (QC-passed reads + QC-failed reads)' \
      '0 + 0 secondary' \
      '0 + 0 supplementary' \
      '0 + 0 duplicates' \
      '900000 + 0 mapped (90.00% : N/A)' \
      '1000000 + 0 paired in sequencing' \
      '500000 + 0 read1' \
      '500000 + 0 read2' \
      '800000 + 0 properly paired (80.00% : N/A)' \
      '850000 + 0 with mate mapped to a different chr' \
      '50000 + 0 with mate mapped to a different chr (mapQ>=5)' \
      > "${prefix}.flagstat"

    touch "${prefix}.stats"
    ${idxstatsStub}
    printf 'id\tbam\texit_code\tstatus\n' > "${prefix}.quickcheck.tsv"
    printf '%s\t%s\t0\tpass\n' "${prefix}" "${bam}" >> "${prefix}.quickcheck.tsv"
    """
}
