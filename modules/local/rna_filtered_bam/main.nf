/*
 * Module: RNA_FILTERED_BAM
 * Runtime command:
 *   bash scripts/core_runtime/RNA_FILTERED_BAM.sh \
 *     <split_name> <solo_dir> <aligned.bam> <canonical_contigs> <outdir> <threads>
 *
 * Inputs:
 *   - STARsolo GeneFull directory from RNA_STARSOLO_ALIGN
 *   - STAR coordinate-sorted aligned BAM from RNA_STARSOLO_ALIGN
 *   - canonical chromosome allowlist resolved once from the STAR index dictionary
 * Outputs:
 *   - low-compression filtered-cells RNA BAM for publication, coverage, and QC
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir } from '../runtime_support/main'

process RNA_FILTERED_BAM {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.rna_filter_retention.tsv"
    publishDir { "${runtimeOutdir()}/rna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*.filtered_cells.bam"

    input:
    tuple val(splitName), val(meta), path(soloDir), path(alignedBam), path(canonicalChromosomes)

    output:
    tuple val(splitName), val(meta), path("${splitName}.filtered_cells.bam"), emit: filtered_bam
    tuple val(splitName), val(meta), path("${splitName}.rna_filter_retention.tsv"), emit: retention_metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock filtered bam for %s\n' "${splitName}" > "${splitName}.filtered_cells.bam"
        pair_count="\$(awk 'NR == 1 { print \$(NF - 1) }' "${alignedBam}")"
        cat > "${splitName}.rna_filter_retention.tsv" <<EOF
split_id	metric	pairs	unit
${splitName}	star_mapped_primary_pairs	\${pair_count}	primary_read1_pair_representatives
${splitName}	paired_filter_pairs	\${pair_count}	primary_read1_pair_representatives
${splitName}	canonical_pairs	\${pair_count}	primary_read1_pair_representatives
${splitName}	called_cell_pairs	\${pair_count}	primary_read1_pair_representatives
EOF

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

        bash "${coreScriptsDir}/RNA_FILTERED_BAM.sh" \\
          "${splitName}" \\
          "${soloDir}" \\
          "${alignedBam}" \\
          "${canonicalChromosomes}" \\
          "." \\
          "${task.cpus}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
