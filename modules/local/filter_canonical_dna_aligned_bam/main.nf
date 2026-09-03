/*
 * Produce the canonical-chromosome aligned BAM used by DNA QC and workflow
 * outputs. GATK MarkDuplicates independently consumes the original alignment,
 * preserving the existing duplicate-marking definition.
 */

process FILTER_CANONICAL_DNA_ALIGNED_BAM {
    tag "${splitName}"
    label 'dna_processing'

    conda "${moduleDir}/../dna_processing/environment-samtools.yml"
    container 'community.wave.seqera.io/library/samtools@sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34'

    input:
    tuple val(splitName), val(meta), path(alignedBam, stageAs: "input.bam"), path(alignedBai, stageAs: "input.bam.bai"), path(canonicalChromosomes, stageAs: "canonical_chromosomes.txt")
    path runtimeScripts, stageAs: 'tresflow/runtime/*'

    output:
    tuple val(splitName), val(meta), path("${splitName}.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}.bam.bai"), emit: bai
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = 'tresflow/runtime'

    if( mode == 'mock' ) {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

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
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        env -u SAMTOOLS_BIN bash "${coreScriptsDir}/FilterCanonicalBam.sh" \\
          "${alignedBam}" \\
          "${splitName}.bam" \\
          "${canonicalChromosomes}" \\
          "${task.cpus}" \\
          normal

        samtools index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}.bam.bai" \\
          "${splitName}.bam"

        chmod a+r "${splitName}.bam" "${splitName}.bam.bai"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
