/*
 * Convert the low-compression filtered RNA BAM into the normally compressed
 * publication artifact. Coverage and QC consume the internal BAM directly on
 * an independent branch.
 */

include { runtimeShellExports; runtimeOutdir } from '../runtime_support/main'

process COMPRESS_RNA_FILTERED_BAM {
    tag "${splitName}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/rna_align" }, mode: 'copy', overwrite: true, pattern: "*.filtered_cells.bam"

    input:
    tuple val(splitName), val(meta), path(internalBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}.filtered_cells.bam"), emit: bam
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def runtimeExports = runtimeShellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        cp -L "${internalBam}" "${splitName}.filtered_cells.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        ${runtimeExports}

        if [[ ! -x "\$SAMTOOLS_BIN" ]]; then
          echo "Missing configured RNA runtime executable: \$SAMTOOLS_BIN" >&2
          exit 1
        fi

        "\$SAMTOOLS_BIN" view \\
          --threads "${task.cpus}" \\
          --bam \\
          --with-header \\
          --output "${splitName}.filtered_cells.bam" \\
          "${internalBam}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
