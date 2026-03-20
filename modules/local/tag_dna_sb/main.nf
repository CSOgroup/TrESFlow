/*
 * Module: TAG_DNA_SAMPLE_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> -D HD=<hd> Tag.codon \
 *     <I2> <R1> <R2> <whitelist> <sample> <tag> <outdir> <first_pass_arg> <rev_comp_arg>
 *
 * Inputs:
 *   - sample metadata
 *   - raw DNA I2 FASTQ as the sample-barcode source
 *   - raw DNA R1 / R2 FASTQs
 *   - shared sample-barcode group map used to derive the effective SB whitelist
 * Outputs:
 *   - DNA sample-barcode-tagged R1 / R2 FASTQs
 *   - barcode counts and summary stats
 */

process TAG_DNA_SAMPLE_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/dna_tagging", mode: 'copy', overwrite: true

    input:
    tuple val(sampleId), val(meta), path(i2), path(r1), path(r2), path(sbGroupMap)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode.R1.fastq"), path("${sampleId}.dna_sample_barcode.R2.fastq"), emit: tagged
    tuple val(sampleId), path("${sampleId}.dna_sample_barcode.counts.tsv"), path("${sampleId}.dna_sample_barcode.stats.tsv"), emit: metrics

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "${params.runtime_python}" "${projectDir}/bin/run_tag.py" \\
      --mode "${mode}" \\
      --script "${params.core_scripts_dir}/Tag.codon" \\
      --i2 "${i2}" \\
      --r1 "${r1}" \\
      --r2 "${r2}" \\
      --sb-group-map "${sbGroupMap}" \\
      --sample "${sampleId}" \\
      --tag "${meta.sample_tag}" \\
      --bc-len ${meta.sample_bc_len} \\
      --bc-start ${meta.sample_bc_start} \\
      --hd ${meta.sample_hd} \\
      --first-pass-arg "${meta.sample_first_pass}" \\
      --rev-comp-arg "${meta.sample_reverse_complement}" \\
      --output-r1 "${sampleId}.dna_sample_barcode.R1.fastq" \\
      --output-r2 "${sampleId}.dna_sample_barcode.R2.fastq" \\
      --output-counts "${sampleId}.dna_sample_barcode.counts.tsv" \\
      --output-stats "${sampleId}.dna_sample_barcode.stats.tsv"
    """
}
