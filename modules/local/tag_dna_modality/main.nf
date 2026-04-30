/*
 * Module: TAG_DNA_MODALITY_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> -D HD=<hd> Tag.codon \
 *     <I2> <tagged_R1> <tagged_R2> <whitelist> <sample> <tag> <outdir> <first_pass_arg> <rev_comp_arg>
 *
 * Inputs:
 *   - sample metadata
 *   - tagmentation-specific DNA index FASTQ as the modality-barcode source
 *   - DNA sample-barcode-tagged R1 / R2 FASTQs
 *   - per-sample DNA modality-barcode whitelist derived from the samplesheet mark mapping
 * Outputs:
 *   - DNA FASTQs tagged with both sample-barcode and modality comments
 *   - modality-barcode counts and summary stats
 */

import RuntimeSupport

process TAG_DNA_MODALITY_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    input:
    tuple val(sampleId), val(meta), path(indexRead), path(taggedR1), path(taggedR2), path(modalityWhitelist)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode_modality.R1.fastq"), path("${sampleId}.dna_sample_barcode_modality.R2.fastq"), emit: tagged
    tuple val(sampleId), path("${sampleId}.dna_modality.counts.tsv"), path("${sampleId}.dna_modality.stats.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')
    def runtimeExports = RuntimeSupport.shellExports(meta)

    """
    ${runtimeExports}

    echo "DNA tagmentation=${meta.dna_tagmentation}; DNA MO index_read=${meta.dna_modality_index_read}; BC_LEN=${meta.modality_bc_len}; BC_START=${meta.modality_bc_start}; HD=${meta.modality_hd}; rev_comp_arg=${meta.modality_reverse_complement}" >&2

    "\$PYTHON3_BIN" "${projectDir}/bin/run_tag.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/Tag.codon" \\
      --i2 "${indexRead}" \\
      --r1 "${taggedR1}" \\
      --r2 "${taggedR2}" \\
      --whitelist "${modalityWhitelist}" \\
      --sample "${sampleId}" \\
      --tag "${meta.modality_tag}" \\
      --bc-len ${meta.modality_bc_len} \\
      --bc-start ${meta.modality_bc_start} \\
      --hd ${meta.modality_hd} \\
      --first-pass-arg "${meta.modality_first_pass}" \\
      --rev-comp-arg "${meta.modality_reverse_complement}" \\
      --output-r1 "${sampleId}.dna_sample_barcode_modality.R1.fastq" \\
      --output-r2 "${sampleId}.dna_sample_barcode_modality.R2.fastq" \\
      --output-counts "${sampleId}.dna_modality.counts.tsv" \\
      --output-stats "${sampleId}.dna_modality.stats.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
