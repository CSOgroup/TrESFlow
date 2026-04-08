/*
 * Module: TAG_RNA_CELL_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D HD=<hd> Tag_Lig3.codon \
 *     <I1> <tagged_R1> <tagged_R2> <whitelist> <sample> <tag> <outdir>
 *
 * Inputs:
 *   - sample metadata
 *   - raw RNA I1 FASTQ as the ligation barcode source
 *   - sample-barcode and UMI-tagged R1 / R2 FASTQs
 *   - cell-barcode whitelist
 * Outputs:
 *   - RNA FASTQs tagged with SB, UM, CB, and RG comments
 *   - per-barcode counts, tag records, and ligation stats
 */

process TAG_RNA_CELL_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    input:
    tuple val(sampleId), val(meta), path(i1), path(taggedR1), path(taggedR2), path(cellWhitelist)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.sample_barcode_umi_cell.R1.fastq.gz"), path("${sampleId}.sample_barcode_umi_cell.R2.fastq.gz"), emit: tagged
    tuple val(sampleId), path("${sampleId}.cell.counts.tsv"), path("${sampleId}.tag_records.tsv"), path("${sampleId}.cell.stats_L1.tsv"), path("${sampleId}.cell.stats_L2.tsv"), path("${sampleId}.cell.stats_L3.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "\$PYTHON3_BIN" "${projectDir}/bin/run_tag_lig3.py" \\
      --mode "${mode}" \\
      --script "${params.core_scripts_dir}/Tag_Lig3.codon" \\
      --i1 "${i1}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --whitelist "${cellWhitelist}" \\
      --sample "${sampleId}" \\
      --tag "${meta.cell_tag}" \\
      --bc-len ${meta.cell_bc_len} \\
      --hd ${meta.cell_hd} \\
      --output-r1 "${sampleId}.sample_barcode_umi_cell.R1.fastq.gz" \\
      --output-r2 "${sampleId}.sample_barcode_umi_cell.R2.fastq.gz" \\
      --output-counts "${sampleId}.cell.counts.tsv" \\
      --output-tag-records "${sampleId}.tag_records.tsv" \\
      --output-stats "${sampleId}.cell.stats_L1.tsv" \\
      --output-stats "${sampleId}.cell.stats_L2.tsv" \\
      --output-stats "${sampleId}.cell.stats_L3.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
