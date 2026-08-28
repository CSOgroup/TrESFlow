/*
 * Module: TAG_DNA_MODALITY_BARCODE
 * Upstream reference:
 *   codon run -plugin seq -release -D BC_LEN=<len> -D BC_START=<start> -D HD=<hd> Tag.codon \
 *     <I2> <tagged_R1> <tagged_R2> <whitelist> <sample> <tag> <outdir> <first_pass_arg> <rev_comp_arg>
 *
 * Inputs:
 *   - sample metadata
 *   - ordered tagmentation-specific DNA index FASTQs as one virtual modality-barcode stream
 *   - DNA sample-barcode-tagged R1 / R2 FASTQs
 *   - per-sample DNA modality-barcode whitelist derived from the samplesheet mark mapping
 * Outputs:
 *   - DNA FASTQs tagged with both sample-barcode and modality comments
 *   - the internal technical-read-set boundary sidecar carried forward unchanged
 *   - modality-barcode counts and summary stats
 */

include { runtimeOutdir } from '../runtime_support/main'

process TAG_DNA_MODALITY_BARCODE {
    tag "${sampleId}"
    label 'codon_wrapper'

    conda "${moduleDir}/../codon_seq/environment.yml"

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.dna_modality.*.tsv"

    input:
    tuple val(sampleId), val(meta), path(indexRead, stageAs: 'index???/*'), path(taggedR1), path(taggedR2), path(modalityWhitelist), path(readSetCounts)
    path helperScripts, stageAs: "tresflow/bin/*"
    path codonScripts, stageAs: "tresflow/codon/*"

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_sample_barcode_modality.R1.fastq"), path("${sampleId}.dna_sample_barcode_modality.R2.fastq"), path(readSetCounts), emit: tagged
    tuple val(sampleId), path("${sampleId}.dna_modality.counts.tsv"), path("${sampleId}.dna_modality.stats.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def indexManifest = (((indexRead instanceof List ? indexRead : [indexRead]).collect { it.toString() }.join('\n')) + '\n').bytes.encodeBase64().toString()
    def taggedR1Manifest = (taggedR1.toString() + '\n').bytes.encodeBase64().toString()
    def taggedR2Manifest = (taggedR2.toString() + '\n').bytes.encodeBase64().toString()

    """
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    printf '%s' '${indexManifest}' | base64 --decode > index.fastq.manifest
    printf '%s' '${taggedR1Manifest}' | base64 --decode > tagged_r1.fastq.manifest
    printf '%s' '${taggedR2Manifest}' | base64 --decode > tagged_r2.fastq.manifest

    echo "DNA tagmentation=${meta.dna_tagmentation}; DNA MO index_read=${meta.dna_modality_index_read}; BC_LEN=${meta.modality_bc_len}; BC_START=${meta.modality_bc_start}; HD=${meta.modality_hd}; rev_comp_arg=${meta.modality_reverse_complement}" >&2

    python3 "tresflow/bin/run_tag.py" \\
      --mode "${mode}" \\
      --script "tresflow/codon/Tag.codon" \\
      --i2-manifest index.fastq.manifest \\
      --r1-manifest tagged_r1.fastq.manifest \\
      --r2-manifest tagged_r2.fastq.manifest \\
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
      --output-stats "${sampleId}.dna_modality.stats.tsv" \\
      --read-set-counts "${readSetCounts}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
