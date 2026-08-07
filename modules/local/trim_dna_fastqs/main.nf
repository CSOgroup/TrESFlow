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

include { runtimeShellExports } from '../runtime_support/main'

process TRIM_DNA_FASTQS {
    tag "${sampleId}"
    label 'codon_wrapper'

    input:
    tuple val(sampleId), val(meta), path(taggedR1), path(taggedR2)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode_modality_cell.R1_val_1.fq"), path("${sampleId}.dna_sample_barcode_modality_cell.R2_val_2.fq"), emit: trimmed
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def runtimeExports = runtimeShellExports(meta)

    """
    ${runtimeExports}

    "\$PYTHON3_BIN" "${projectDir}/bin/run_trim_galore.py" \\
      --mode "${mode}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --trim-galore-bin "\$TRIM_GALORE_BIN" \\
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
