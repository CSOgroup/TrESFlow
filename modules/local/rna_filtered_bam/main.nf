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

include { runtimeOutdir } from '../runtime_support/main'

process RNA_FILTERED_BAM {
    tag "${splitName}"
    label 'rna_alignment'

    conda "${moduleDir}/../rna_alignment/environment-filter.yml"
    container 'community.wave.seqera.io/library/samtools_python@sha256:74535d380b6c327aa8a82ad941f00d900d26f1a74217e82679f9d64b1b9e28d3'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.rna_filter_retention.tsv"
    publishDir { "${runtimeOutdir()}/rna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*.filtered_cells.bam"

    input:
    tuple val(splitName), val(meta), path(soloDir), path(alignedBam), path(canonicalChromosomes)
    path runtimeScripts, stageAs: 'tresflow/runtime/*'

    output:
    tuple val(splitName), val(meta), path("${splitName}.filtered_cells.bam"), emit: filtered_bam
    tuple val(splitName), val(meta), path("${splitName}.rna_filter_retention.tsv"), emit: retention_metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = 'tresflow/runtime'
    if( mode == 'mock' ) {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

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
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        bash "${coreScriptsDir}/RNA_FILTERED_BAM.sh" \\
          "${splitName}" \\
          "${soloDir}" \\
          "${alignedBam}" \\
          "${canonicalChromosomes}" \\
          "." \\
          "${task.cpus}"

        chmod a+r "${splitName}.filtered_cells.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
