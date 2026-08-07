/*
 * Produce the canonical-chromosome aligned BAM used by DNA QC and workflow
 * outputs. GATK MarkDuplicates independently consumes the original alignment,
 * preserving the existing duplicate-marking definition.
 */

include { runtimeShellExports; runtimeCoreScriptsDir } from '../runtime_support/main'

process FILTER_CANONICAL_DNA_ALIGNED_BAM {
    tag "${splitName}"
    label 'process_single'

    input:
    tuple val(splitName), val(meta), path(alignedBam, stageAs: "input.bam"), path(alignedBai, stageAs: "input.bam.bai"), path(canonicalChromosomes, stageAs: "canonical_chromosomes.txt")

    output:
    tuple val(splitName), val(meta), path("${splitName}.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}.bam.bai"), emit: bai
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        cp -L "${alignedBam}" "${splitName}.bam"
        cp -L "${alignedBai}" "${splitName}.bam.bai"

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
          "${alignedBam}" \\
          "${splitName}.bam" \\
          "${canonicalChromosomes}" \\
          "${task.cpus}" \\
          normal

        "\$SAMTOOLS_BIN" index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}.bam.bai" \\
          "${splitName}.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
