/*
 * Module: TAG_DNA_CELL_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D HD=<hd> Tag_Lig3.codon \
 *     <I1> <tagged_R1> <tagged_R2> <whitelist> <sample> <tag> <outdir>
 *
 * Inputs:
 *   - sample metadata
 *   - raw DNA I1 FASTQ as the ligation barcode source
 *   - DNA sample-barcode and modality-tagged R1 / R2 FASTQs
 *   - DNA ligation whitelist
 * Outputs:
 *   - DNA FASTQs tagged with SB, MO, CB, and RG comments
 *   - per-barcode counts, tag records, and ligation stats
 */

process TAG_DNA_CELL_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/dna_tagging", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId), val(meta), path(i1), path(taggedR1), path(taggedR2), path(cellWhitelist)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode_modality_cell.R1.fastq"), path("${sampleId}.dna_sample_barcode_modality_cell.R2.fastq"), emit: tagged
    tuple val(sampleId), path("${sampleId}.dna_cell.counts.tsv"), path("${sampleId}.dna_tag_records.tsv"), path("${sampleId}.dna_cell.stats_L1.tsv"), path("${sampleId}.dna_cell.stats_L2.tsv"), path("${sampleId}.dna_cell.stats_L3.tsv"), emit: metrics

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "${params.runtime_python}" "${projectDir}/bin/run_tag_lig3.py" \\
      --mode "${mode}" \\
      --script "${params.upstream_dir}/Tag_Lig3.codon" \\
      --i1 "${i1}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --whitelist "${cellWhitelist}" \\
      --sample "${sampleId}" \\
      --tag "${meta.cell_tag}" \\
      --bc-len ${meta.cell_bc_len} \\
      --hd ${meta.cell_hd} \\
      --output-r1 "${sampleId}.dna_sample_barcode_modality_cell.R1.fastq" \\
      --output-r2 "${sampleId}.dna_sample_barcode_modality_cell.R2.fastq" \\
      --output-counts "${sampleId}.dna_cell.counts.tsv" \\
      --output-tag-records "${sampleId}.dna_tag_records.tsv" \\
      --output-stats "${sampleId}.dna_cell.stats_L1.tsv" \\
      --output-stats "${sampleId}.dna_cell.stats_L2.tsv" \\
      --output-stats "${sampleId}.dna_cell.stats_L3.tsv"
    """
}
