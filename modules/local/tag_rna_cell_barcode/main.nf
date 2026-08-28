/*
 * Module: TAG_RNA_CELL_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D HD=<hd> Tag_Lig3.codon \
 *     <I1> <tagged_R1> <tagged_R2> <whitelist> <sample> <tag> <outdir>
 *
 * Inputs:
 *   - sample metadata
 *   - ordered raw RNA I1 FASTQs as one virtual ligation-barcode stream
 *   - sample-barcode and UMI-tagged R1 / R2 FASTQs
 *   - internal technical-read-set boundary counts from the first tag stage
 *   - cell-barcode whitelist
 * Outputs:
 *   - RNA FASTQs tagged with SB, UM, CB, and RG comments
 *   - per-barcode counts, tag records, and ligation stats
 */

include { runtimeOutdir } from '../runtime_support/main'

process TAG_RNA_CELL_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    conda "${moduleDir}/../codon_seq/environment.yml"

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.rna_cell.*.tsv"
    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.rna_tag_records.tsv.gz"

    input:
    tuple val(sampleId), val(meta), path(i1, stageAs: 'i1???/*'), path(taggedR1), path(taggedR2), path(cellWhitelist), path(readSetCounts)
    path helperScripts, stageAs: "tresflow/bin/*"
    path codonScripts, stageAs: "tresflow/codon/*"

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.sample_barcode_umi_cell.R1.fastq"), path("${sampleId}.sample_barcode_umi_cell.R2.fastq"), emit: tagged
    tuple val(sampleId), path("${sampleId}.cell.counts.tsv"), path("${sampleId}.cell.stats_L1.tsv"), path("${sampleId}.cell.stats_L2.tsv"), path("${sampleId}.cell.stats_L3.tsv"), emit: metrics
    path("${sampleId}.rna_cell.*.tsv"), emit: tres_cell_stats
    path("${sampleId}.rna_tag_records.tsv.gz"), emit: tres_tag_records
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def i1Manifest = (((i1 instanceof List ? i1 : [i1]).collect { it.toString() }.join('\n')) + '\n').bytes.encodeBase64().toString()
    def taggedR1Manifest = (taggedR1.toString() + '\n').bytes.encodeBase64().toString()
    def taggedR2Manifest = (taggedR2.toString() + '\n').bytes.encodeBase64().toString()

    """
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    printf '%s' '${i1Manifest}' | base64 --decode > i1.fastq.manifest
    printf '%s' '${taggedR1Manifest}' | base64 --decode > tagged_r1.fastq.manifest
    printf '%s' '${taggedR2Manifest}' | base64 --decode > tagged_r2.fastq.manifest

    python3 "tresflow/bin/run_tag_lig3.py" \\
      --mode "${mode}" \\
      --script "tresflow/codon/Tag_Lig3.codon" \\
      --i1-manifest i1.fastq.manifest \\
      --r1-manifest tagged_r1.fastq.manifest \\
      --r2-manifest tagged_r2.fastq.manifest \\
      --whitelist "${cellWhitelist}" \\
      --sample "${sampleId}" \\
      --tag "${meta.cell_tag}" \\
      --bc-len ${meta.cell_bc_len} \\
      --hd ${meta.cell_hd} \\
      --output-r1 "${sampleId}.sample_barcode_umi_cell.R1.fastq" \\
      --output-r2 "${sampleId}.sample_barcode_umi_cell.R2.fastq" \\
      --output-counts "${sampleId}.cell.counts.tsv" \\
      --output-tag-records "${sampleId}.rna_tag_records.tsv" \\
      --read-set-counts "${readSetCounts}" \\
      --output-stats "${sampleId}.cell.stats_L1.tsv" \\
      --output-stats "${sampleId}.cell.stats_L2.tsv" \\
      --output-stats "${sampleId}.cell.stats_L3.tsv"

    cp "${sampleId}.cell.counts.tsv" "${sampleId}.rna_cell.counts.tsv"
    cp "${sampleId}.cell.stats_L1.tsv" "${sampleId}.rna_cell.stats_L1.tsv"
    cp "${sampleId}.cell.stats_L2.tsv" "${sampleId}.rna_cell.stats_L2.tsv"
    cp "${sampleId}.cell.stats_L3.tsv" "${sampleId}.rna_cell.stats_L3.tsv"

    pigz -f -p "${task.cpus}" "${sampleId}.rna_tag_records.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
