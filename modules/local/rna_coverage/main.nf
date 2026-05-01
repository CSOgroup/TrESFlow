/*
 * Module: RNA_COVERAGE
 * Runtime command:
 *   bash scripts/core_runtime/RNA_COVERAGE.sh \
 *     <split_name> <filtered_cells.bam> <star_index_dir> <chrom.sizes> <outdir> <threads>
 *
 * Inputs:
 *   - filtered-cells RNA BAM from RNA_FILTERED_BAM
 *   - exact RNA STAR index directory from references.rna_ref_dir and derived STAR chrNameLength.txt chromosome sizes file
 * Outputs:
 *   - stranded and unstranded RNA bigWig tracks
 */

import RuntimeSupport

process RNA_COVERAGE {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/rna_align", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(filteredBam), val(starIndexDir), val(chromSizes)

    output:
    tuple val(splitName), val(meta), path("${splitName}.stranded_*.bw"), optional: true, emit: stranded_bw
    tuple val(splitName), val(meta), path("${splitName}.unstranded_*.bw"), optional: true, emit: unstranded_bw
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')
    def runtimeExports = RuntimeSupport.shellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock stranded bigwig\n' > "${splitName}.stranded_Signal.Unique.str1.out.bw"
        printf 'mock unstranded bigwig\n' > "${splitName}.unstranded_Signal.Unique.str1.out.bw"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        ${runtimeExports}

        for required_bin in "\$STAR_BIN" "\$BEDGRAPH_TO_BIGWIG_BIN"; do
          if [[ ! -x "\${required_bin}" ]]; then
            echo "Missing configured RNA runtime executable: \${required_bin}" >&2
            exit 1
          fi
        done

        bash "${coreScriptsDir}/RNA_COVERAGE.sh" \\
          "${splitName}" \\
          "${filteredBam}" \\
          "${starIndexDir}" \\
          "${chromSizes}" \\
          "." \\
          "${task.cpus}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
