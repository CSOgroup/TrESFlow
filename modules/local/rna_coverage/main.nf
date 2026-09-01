/*
 * Module: RNA_COVERAGE
 * Runtime command:
 *   bash scripts/core_runtime/RNA_COVERAGE.sh \
 *     <split_name> <filtered_cells.bam> <star_index_dir> <chrom.sizes> <outdir> <threads>
 *
 * Inputs:
 *   - filtered-cells RNA BAM from RNA_FILTERED_BAM
 *   - exact RNA STAR index directory and canonical-only chromosome sizes resolved from its dictionary
 * Outputs:
 *   - stranded and unstranded RNA bigWig tracks
 */

include { runtimeOutdir } from '../runtime_support/main'

process RNA_COVERAGE {
    tag "${splitName}"
    label 'rna_alignment'

    conda "${moduleDir}/../rna_alignment/environment-coverage.yml"
    container 'community.wave.seqera.io/library/star_ucsc-bedgraphtobigwig@sha256:133ac55ecc30285f3e8b0efa7e3577efd6184354641098d0a290c66268038d88'

    publishDir { "${runtimeOutdir()}/rna_align" }, mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(filteredBam), path(starIndexDir), path(chromSizes)
    path runtimeScripts, stageAs: 'tresflow/runtime/*'

    output:
    tuple val(splitName), val(meta), path("${splitName}.stranded_*.bw"), optional: true, emit: stranded_bw
    tuple val(splitName), val(meta), path("${splitName}.unstranded_*.bw"), optional: true, emit: unstranded_bw
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = 'tresflow/runtime'
    if( mode == 'mock' ) {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        printf 'mock stranded bigwig\n' > "${splitName}.stranded_Signal.Unique.str1.out.bw"
        printf 'mock unstranded bigwig\n' > "${splitName}.unstranded_Signal.Unique.str1.out.bw"

        printf '        "%s":\n          component: "local"\n        END_VERSIONS\n' \
          "${task.process}" > versions.yml
        """
    }
    else {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        bash "${coreScriptsDir}/RNA_COVERAGE.sh" \\
          "${splitName}" \\
          "${filteredBam}" \\
          "${starIndexDir}" \\
          "${chromSizes}" \\
          "." \\
          "${task.cpus}"

        printf '        "%s":\n          component: "local"\n        END_VERSIONS\n' \
          "${task.process}" > versions.yml
        """
    }
}
