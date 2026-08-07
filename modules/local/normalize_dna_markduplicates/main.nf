/*
 * Normalize nf-core/gatk4/markduplicates outputs back to the TrESFlow DNA
 * filename contract and remove noncanonical alignments after duplicate marking.
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir } from '../runtime_support/main'

process NORMALIZE_DNA_MARKDUPLICATES {
    tag "${splitName}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/dna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*_MarkedDup.bam*"
    publishDir { "${runtimeOutdir()}/dna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*.DuplicateMetrics.txt"

    input:
    tuple val(splitName), val(meta), path(markedDupBam, stageAs: "input_markeddup.bam"), path(markedDupBai, stageAs: "input_markeddup.bam.bai"), path(markedDupMetrics, stageAs: "input.DuplicateMetrics.txt"), path(canonicalChromosomes, stageAs: "canonical_chromosomes.txt")

    output:
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}.DuplicateMetrics.txt"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

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
        ${runtimeExports}

        bash "${coreScriptsDir}/FilterCanonicalBam.sh" \\
          "${markedDupBam}" \\
          "${splitName}_MarkedDup.bam" \\
          "${canonicalChromosomes}" \\
          "${task.cpus}" \\
          normal

        "\$SAMTOOLS_BIN" index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}_MarkedDup.bam.bai" \\
          "${splitName}_MarkedDup.bam"

        cp -L "${markedDupMetrics}" "${splitName}.DuplicateMetrics.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
