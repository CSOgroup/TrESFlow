/*
 * Module: TAG_RNA_SAMPLE_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> -D HD=<hd> Tag.codon \
 *     <I2> <R1> <R2> <whitelist> <sample> <tag> <outdir> <first_pass_arg> <rev_comp_arg>
 *
 * This wrapper keeps the upstream command as the default path and exposes stable output names
 * for Nextflow integration. In the test profile it switches to a lightweight mock implementation
 * so `-profile test` runs without Codon.
 *
 * Inputs:
 *   - sample metadata
 *   - raw RNA R1 FASTQ
 *   - raw RNA R2 FASTQ (used both as read 2 and as the barcode source for this RNA step)
 *   - RNA SB-group map TSV used as the single source of truth for experiment-used sample barcodes
 * Outputs:
 *   - sample-barcode-tagged R1 / R2 FASTQs
 *   - barcode counts and summary stats
 */

process TAG_RNA_SAMPLE_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    input:
    tuple val(sampleId), val(meta), path(r1), path(r2), path(sbGroupMap)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.sample_barcode.R1.fastq.gz"), path("${sampleId}.sample_barcode.R2.fastq.gz"), emit: tagged
    tuple val(sampleId), path("${sampleId}.sample_barcode.counts.tsv"), path("${sampleId}.sample_barcode.stats.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    "\$PYTHON3_BIN" "${projectDir}/bin/run_tag.py" \\
      --mode "${mode}" \\
      --script "${params.core_scripts_dir}/Tag.codon" \\
      --i2 "${r2}" \\
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
      --output-r1 "${sampleId}.sample_barcode.R1.fastq.gz" \\
      --output-r2 "${sampleId}.sample_barcode.R2.fastq.gz" \\
      --output-counts "${sampleId}.sample_barcode.counts.tsv" \\
      --output-stats "${sampleId}.sample_barcode.stats.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
