/*
 * Module: TAG_RNA_UMI
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> Tag_UMI.codon \
 *     <I2> <tagged_R1> <tagged_R2> <sample> <tag> <outdir>
 *
 * Inputs:
 *   - sample metadata
 *   - ordered raw RNA R2 FASTQs as one virtual UMI-barcode stream
 *   - sample-barcode-tagged R1 / R2 FASTQs from TAG_RNA_SAMPLE_BARCODE
 * Outputs:
 *   - RNA FASTQs tagged with both sample-barcode and UMI comments
 *   - the internal technical-read-set boundary sidecar carried forward unchanged
 *   - UMI counts table
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir } from '../runtime_support/main'

process TAG_RNA_UMI {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.rna_umi.counts.tsv"

    input:
    tuple val(sampleId), val(meta), path(rawR2, stageAs: 'raw_r2???/*'), path(taggedR1), path(taggedR2), path(readSetCounts)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.sample_barcode_umi.R1.fastq"), path("${sampleId}.sample_barcode_umi.R2.fastq"), path(readSetCounts), emit: tagged
    tuple val(sampleId), path("${sampleId}.umi.counts.tsv"), emit: metrics
    path("${sampleId}.rna_umi.counts.tsv"), emit: tres_stats
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)
    def rawR2Manifest = (((rawR2 instanceof List ? rawR2 : [rawR2]).collect { it.toString() }.join('\n')) + '\n').bytes.encodeBase64().toString()
    def taggedR1Manifest = (taggedR1.toString() + '\n').bytes.encodeBase64().toString()
    def taggedR2Manifest = (taggedR2.toString() + '\n').bytes.encodeBase64().toString()

    """
    ${runtimeExports}

    printf '%s' '${rawR2Manifest}' | base64 --decode > raw_r2.fastq.manifest
    printf '%s' '${taggedR1Manifest}' | base64 --decode > tagged_r1.fastq.manifest
    printf '%s' '${taggedR2Manifest}' | base64 --decode > tagged_r2.fastq.manifest

    "\$PYTHON3_BIN" "${projectDir}/bin/run_tag_umi.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/Tag_UMI.codon" \\
      --i2-manifest raw_r2.fastq.manifest \\
      --r1-manifest tagged_r1.fastq.manifest \\
      --r2-manifest tagged_r2.fastq.manifest \\
      --sample "${sampleId}" \\
      --tag "${meta.umi_tag}" \\
      --bc-len ${meta.umi_bc_len} \\
      --bc-start ${meta.umi_bc_start} \\
      --output-r1 "${sampleId}.sample_barcode_umi.R1.fastq" \\
      --output-r2 "${sampleId}.sample_barcode_umi.R2.fastq" \\
      --output-counts "${sampleId}.umi.counts.tsv" \\
      --read-set-counts "${readSetCounts}" \\
      --rev-comp

    cp "${sampleId}.umi.counts.tsv" "${sampleId}.rna_umi.counts.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
