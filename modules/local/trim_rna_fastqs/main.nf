/*
 * Module: TRIM_RNA_FASTQS
 * Upstream reference:
 *   trim_galore \
 *     --quality 10 \
 *     --cores <task.cpus> \
 *     --output_dir <outdir> \
 *     --gzip \
 *     --length 20 \
 *     --paired \
 *     <CB_tagged_R1.fastq> <CB_tagged_R2.fastq>
 *
 * Inputs:
 *   - sample metadata
 *   - RNA R1 / R2 FASTQs tagged with SB, UM, CB, and RG comments
 * Outputs:
 *   - trim_galore paired-end FASTQs named with the standard _val_1 / _val_2 suffixes
 *
 * Notes:
 *   - The upstream launcher deletes the untrimmed CB FASTQs after trimming.
 *     This Nextflow slice keeps those earlier module outputs published and advances on the trimmed files.
 */

process TRIM_RNA_FASTQS {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/tagging", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId), val(meta), path(taggedR1), path(taggedR2)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.sample_barcode_umi_cell.R1_val_1.fq.gz"), path("${sampleId}.sample_barcode_umi_cell.R2_val_2.fq.gz"), emit: trimmed

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "\$PYTHON3_BIN" "${projectDir}/bin/run_trim_galore.py" \\
      --mode "${mode}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --trim-galore-bin "\$TRIM_GALORE_BIN" \\
      --quality 10 \\
      --cores ${task.cpus} \\
      --length 20 \\
      --output-r1 "${sampleId}.sample_barcode_umi_cell.R1_val_1.fq.gz" \\
      --output-r2 "${sampleId}.sample_barcode_umi_cell.R2_val_2.fq.gz"
    """
}
