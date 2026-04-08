/*
 * Module: SPLIT_RNA_READS
 * Upstream reference:
 *   codon run -plugin seq -release Split_ReadsV2.codon \
 *     <Sample> <OutFolder> <LibName> rna - <trimmed_R1.fq.gz> <trimmed_R2.fq.gz> <sb_group_map.tsv>
 *
 * Inputs:
 *   - sample metadata
 *   - trim_galore RNA FASTQs from the CB-tagged reads
 *   - shared sample-barcode group map TSV keyed by sample and group
 * Outputs:
 *   - per-group RNA FASTQ pairs named as upstream Split_ReadsV2 outputs
 *   - per-group SAM RG header TSVs named as upstream Split_ReadsV2 outputs
 *
 * Notes:
 *   - The upstream sample-barcode group map example uses full SB strings even though the script comments
 *     discuss dropping an injected leading base. This wrapper follows the actual script logic:
 *     raw SB match first, then drop-first fallback.
 */

process SPLIT_RNA_READS {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/split", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId), val(meta), path(trimmedR1), path(trimmedR2), path(sbGroupMap)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fq.gz"), path("${sampleId}_*_R2.fq.gz"), emit: split_fastqs
    tuple val(sampleId), val(meta), path("SAM_RG_Header_${sampleId}_*.tsv"), emit: rg_headers
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "\$PYTHON3_BIN" "${projectDir}/bin/run_split_reads_rna.py" \\
      --mode "${mode}" \\
      --script "${params.core_scripts_dir}/Split_ReadsV2.codon" \\
      --r1 "${trimmedR1}" \\
      --r2 "${trimmedR2}" \\
      --sb-group-map "${sbGroupMap}" \\
      --sample "${sampleId}" \\
      --library-name "${meta.library_name}" \\
      --output-dir "."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
