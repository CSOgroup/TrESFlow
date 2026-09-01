/*
 * Module: TRIM_DNA_FASTQS
 * Upstream reference:
 *   trim_galore \
 *     --quality 10 \
 *     --cores <task.cpus> \
 *     --output_dir <outdir> \
 *     --dont_gzip \
 *     --length 20 \
 *     --paired \
 *     <CB_tagged_R1.fastq> <CB_tagged_R2.fastq>
 *
 * Inputs:
 *   - sample metadata
 *   - DNA R1 / R2 FASTQs tagged with SB, MO, CB, and RG comments
 * Outputs:
 *   - uncompressed trim_galore paired-end FASTQs named with the standard _val_1 / _val_2 suffixes
 */

process TRIM_DNA_FASTQS {
    tag "${sampleId}"
    label 'fastq_preprocessing'

    conda "${moduleDir}/../fastq_preprocessing/environment-trim.yml"
    container 'quay.io/biocontainers/trim-galore@sha256:a02bb87b8ce02d86efd0ffd65e2cce1559b52689faab42faad1df145657390cf'

    input:
    tuple val(sampleId), val(meta), path(taggedR1), path(taggedR2)
    path helperScript, stageAs: 'tresflow/bin/run_trim_galore.py'

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode_modality_cell.R1_val_1.fq"), path("${sampleId}.dna_sample_barcode_modality_cell.R2_val_2.fq"), emit: trimmed
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    python3 "tresflow/bin/run_trim_galore.py" \\
      --mode "${mode}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --quality 10 \\
      --cores ${task.cpus} \\
      --length 20 \\
      --output-r1 "${sampleId}.dna_sample_barcode_modality_cell.R1_val_1.fq" \\
      --output-r2 "${sampleId}.dna_sample_barcode_modality_cell.R2_val_2.fq"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
