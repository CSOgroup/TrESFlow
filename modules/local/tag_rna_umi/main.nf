/*
 * Module: TAG_RNA_UMI
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> Tag_UMI.codon \
 *     <I2> <tagged_R1> <tagged_R2> <sample> <tag> <outdir>
 *
 * Inputs:
 *   - sample metadata
 *   - raw RNA R2 FASTQ as the UMI barcode source
 *   - sample-barcode-tagged R1 / R2 FASTQs from TAG_RNA_SAMPLE_BARCODE
 * Outputs:
 *   - RNA FASTQs tagged with both sample-barcode and UMI comments
 *   - UMI counts table
 */

process TAG_RNA_UMI {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/tagging", mode: 'copy', overwrite: true, saveAs: { filename ->
        filename ==~ /.*\.(fastq|fq)(\.gz)?$/ ? null : filename
    }

    input:
    tuple val(sampleId), val(meta), path(rawR2), path(taggedR1), path(taggedR2)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.sample_barcode_umi.R1.fastq.gz"), path("${sampleId}.sample_barcode_umi.R2.fastq.gz"), emit: tagged
    tuple val(sampleId), path("${sampleId}.umi.counts.tsv"), emit: metrics

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "\$PYTHON3_BIN" "${projectDir}/bin/run_tag_umi.py" \\
      --mode "${mode}" \\
      --script "${params.core_scripts_dir}/Tag_UMI.codon" \\
      --i2 "${rawR2}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --sample "${sampleId}" \\
      --tag "${meta.umi_tag}" \\
      --bc-len ${meta.umi_bc_len} \\
      --bc-start ${meta.umi_bc_start} \\
      --output-r1 "${sampleId}.sample_barcode_umi.R1.fastq.gz" \\
      --output-r2 "${sampleId}.sample_barcode_umi.R2.fastq.gz" \\
      --output-counts "${sampleId}.umi.counts.tsv" \\
      --rev-comp
    """
}
