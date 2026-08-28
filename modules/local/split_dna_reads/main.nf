/*
 * Module: SPLIT_DNA_READS
 * Upstream reference:
 *   codon run -plugin seq -release Split_ReadsV2.codon \
 *     <Sample> <OutFolder> <LibName> dna <mo_map.tsv> <split_input_R1.fq> <split_input_R2.fq> <sb_group_map.tsv>
 *
 * Inputs:
 *   - sample metadata
 *   - uncompressed paired DNA FASTQs entering splitting: Trim Galore outputs
 *     directly, or post-artifact-filter outputs for enabled dual-tag samples
 *   - DNA modality map TSV keyed by sample, group, mark, and modality barcode
 *   - shared sample-barcode group map TSV keyed by sample and group
 * Outputs:
 *   - uncompressed per-group per-mark DNA FASTQ pairs for downstream computation
 *   - per-group per-mark SAM RG header TSVs named as upstream Split_ReadsV2 outputs
 */

include { runtimeOutdir } from '../runtime_support/main'

process SPLIT_DNA_READS {
    tag "${sampleId}"
    label 'codon_wrapper'

    conda "${moduleDir}/../codon_seq/environment.yml"

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.dna_read_retention.tsv"

    input:
    tuple val(sampleId), val(meta), path(splitInputR1), path(splitInputR2), path(moMap), path(sbGroupMap)
    path helperScripts, stageAs: "tresflow/bin/*"
    path codonScripts, stageAs: "tresflow/codon/*"

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fastq"), path("${sampleId}_*_R2.fastq"), emit: split_fastqs
    tuple val(sampleId), val(meta), path("SAM_RG_Header_${sampleId}_*.tsv"), emit: rg_headers
    tuple val(sampleId), val(meta), path("${sampleId}.dna_read_retention.tsv"), emit: retention_metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    """
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    python3 "tresflow/bin/run_split_reads_dna.py" \\
      --mode "${mode}" \\
      --script "tresflow/codon/Split_ReadsV2.codon" \\
      --r1 "${splitInputR1}" \\
      --r2 "${splitInputR2}" \\
      --mo-map "${moMap}" \\
      --sb-group-map "${sbGroupMap}" \\
      --sample "${sampleId}" \\
      --library-name "${meta.library_name}" \\
      --output-dir "."

    printf '%s\\n' \\
      '"${task.process}":' \\
      '  component: "local"' \\
      > versions.yml
    """
}
