/*
 * Module: TAG_DNA_SAMPLE_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> -D HD=<hd> Tag.codon \
 *     <I2> <R1> <R2> <whitelist> <sample> <tag> <outdir> <first_pass_arg> <rev_comp_arg>
 *
 * Inputs:
 *   - sample metadata
 *   - tagmentation-specific DNA index FASTQ as the sample-barcode source
 *   - raw DNA R1 / R2 FASTQs
 *   - shared sample-barcode group map used to derive the effective SB whitelist
 * Outputs:
 *   - DNA sample-barcode-tagged R1 / R2 FASTQs
 *   - barcode counts and summary stats
 */

import RuntimeSupport

process TAG_DNA_SAMPLE_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    input:
    tuple val(sampleId), val(meta), path(indexRead), path(r1), path(r2), path(sbGroupMap)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode.R1.fastq"), path("${sampleId}.dna_sample_barcode.R2.fastq"), emit: tagged
    tuple val(sampleId), path("${sampleId}.dna_sample_barcode.counts.tsv"), path("${sampleId}.dna_sample_barcode.stats.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')
    def runtimeExports = RuntimeSupport.shellExports(meta)

    """
    ${runtimeExports}

    echo "DNA tagmentation=${meta.dna_tagmentation}; DNA SB source=${meta.dna_sb_barcode_source}; DNA SB length=${meta.dna_sb_barcode_len}; index_read=${meta.dna_sample_index_read}; BC_LEN=${meta.sample_bc_len}; BC_START=${meta.sample_bc_start}; HD=${meta.sample_hd}; rev_comp_arg=${meta.sample_reverse_complement}" >&2

    "\$PYTHON3_BIN" "${projectDir}/bin/run_tag.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/Tag.codon" \\
      --i2 "${indexRead}" \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
