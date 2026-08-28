/*
 * Module: TAG_DNA_SAMPLE_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> -D HD=<hd> Tag.codon \
 *     <I2> <R1> <R2> <whitelist> <sample> <tag> <outdir> <first_pass_arg> <rev_comp_arg>
 *
 * Inputs:
 *   - sample metadata
 *   - ordered tagmentation-specific DNA index FASTQs as one virtual sample-barcode stream
 *   - ordered raw DNA R1 / R2 FASTQ collections
 *   - shared sample-barcode group map used to derive the effective SB whitelist
 * Outputs:
 *   - DNA sample-barcode-tagged R1 / R2 FASTQs
 *   - internal technical-read-set record counts used to preserve boundaries downstream
 *   - barcode counts and summary stats
 */

include { runtimeOutdir } from '../runtime_support/main'

process TAG_DNA_SAMPLE_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    conda "${moduleDir}/../codon_seq/environment.yml"

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.dna_sample_barcode.*.tsv"

    input:
    tuple val(sampleId), val(meta), path(indexRead, stageAs: 'index???/*'), path(r1, stageAs: 'r1???/*'), path(r2, stageAs: 'r2???/*'), path(sbGroupMap)
    path helperScripts, stageAs: "tresflow/bin/*"
    path codonScripts, stageAs: "tresflow/codon/*"

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode.R1.fastq"), path("${sampleId}.dna_sample_barcode.R2.fastq"), path("${sampleId}.dna.read_set_counts.tsv"), emit: tagged
    tuple val(sampleId), path("${sampleId}.dna_sample_barcode.counts.tsv"), path("${sampleId}.dna_sample_barcode.stats.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def indexManifest = (((indexRead instanceof List ? indexRead : [indexRead]).collect { it.toString() }.join('\n')) + '\n').bytes.encodeBase64().toString()
    def r1Manifest = (((r1 instanceof List ? r1 : [r1]).collect { it.toString() }.join('\n')) + '\n').bytes.encodeBase64().toString()
    def r2Manifest = (((r2 instanceof List ? r2 : [r2]).collect { it.toString() }.join('\n')) + '\n').bytes.encodeBase64().toString()

    """
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    printf '%s' '${indexManifest}' | base64 --decode > index.fastq.manifest
    printf '%s' '${r1Manifest}' | base64 --decode > r1.fastq.manifest
    printf '%s' '${r2Manifest}' | base64 --decode > r2.fastq.manifest

    echo "DNA tagmentation=${meta.dna_tagmentation}; DNA SB source=${meta.dna_sb_barcode_source}; DNA SB length=${meta.dna_sb_barcode_len}; index_read=${meta.dna_sample_index_read}; BC_LEN=${meta.sample_bc_len}; BC_START=${meta.sample_bc_start}; HD=${meta.sample_hd}; rev_comp_arg=${meta.sample_reverse_complement}" >&2

    python3 "tresflow/bin/run_tag.py" \\
      --mode "${mode}" \\
      --script "tresflow/codon/Tag.codon" \\
      --i2-manifest index.fastq.manifest \\
      --r1-manifest r1.fastq.manifest \\
      --r2-manifest r2.fastq.manifest \\
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
      --output-stats "${sampleId}.dna_sample_barcode.stats.tsv" \\
      --output-read-set-counts "${sampleId}.dna.read_set_counts.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
