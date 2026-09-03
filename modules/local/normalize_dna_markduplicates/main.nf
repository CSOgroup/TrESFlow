/*
 * Normalize nf-core/gatk4/markduplicates outputs back to the TrESFlow DNA
 * filename contract and remove noncanonical alignments after duplicate marking.
 */

include { runtimeOutdir } from '../runtime_support/main'

process NORMALIZE_DNA_MARKDUPLICATES {
    tag "${splitName}"
    label 'dna_processing'

    conda "${moduleDir}/../dna_processing/environment-samtools.yml"
    container 'community.wave.seqera.io/library/samtools@sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34'

    publishDir { "${runtimeOutdir()}/dna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*_MarkedDup.bam*"
    publishDir { "${runtimeOutdir()}/dna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*.DuplicateMetrics.txt"

    input:
    tuple val(splitName), val(meta), path(markedDupBam, stageAs: "input_markeddup.bam"), path(markedDupBai, stageAs: "input_markeddup.bam.bai"), path(markedDupMetrics, stageAs: "input.DuplicateMetrics.txt"), path(canonicalChromosomes, stageAs: "canonical_chromosomes.txt")
    path runtimeScripts, stageAs: 'tresflow/runtime/*'

    output:
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}.DuplicateMetrics.txt"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = 'tresflow/runtime'

    if( mode == 'mock' ) {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        cp -L "${markedDupBam}" "${splitName}_MarkedDup.bam"
        cp -L "${markedDupBai}" "${splitName}_MarkedDup.bam.bai"
        cp -L "${markedDupMetrics}" "${splitName}.DuplicateMetrics.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        env -u SAMTOOLS_BIN bash "${coreScriptsDir}/FilterCanonicalBam.sh" \\
          "${markedDupBam}" \\
          "${splitName}_MarkedDup.bam" \\
          "${canonicalChromosomes}" \\
          "${task.cpus}" \\
          normal

        samtools index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}_MarkedDup.bam.bai" \\
          "${splitName}_MarkedDup.bam"

        chmod a+r "${splitName}_MarkedDup.bam" "${splitName}_MarkedDup.bam.bai"

        cp -L "${markedDupMetrics}" "${splitName}.DuplicateMetrics.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
