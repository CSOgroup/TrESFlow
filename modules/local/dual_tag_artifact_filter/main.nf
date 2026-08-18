/*
 * Module: DUAL_TAG_ARTIFACT_FILTER
 *
 * Discard dual-tagmentation DNA read pairs when either mate contains one of
 * the audited 48 exact 23-nt linker signatures remaining after Trim Galore.
 * Cutadapt uses action=none so retained trimmed FASTQ records are passed
 * through without sequence modification. This process is routed only
 * dual-tagmentation samples by DNA_CORE.
 */

include { runtimeShellExports; runtimeOutdir } from '../runtime_support/main'

process DUAL_TAG_ARTIFACT_FILTER {
    tag "${sampleId}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.dual_tag_artifact_filter.cutadapt.json"
    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.dual_tag_artifact_filter.summary.tsv"

    input:
    tuple val(sampleId), val(meta), path(trimmedR1), path(trimmedR2)
    path signatureFasta

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.dna_dual_tag_clean.R1.fastq"), path("${sampleId}.dna_dual_tag_clean.R2.fastq"), emit: filtered
    tuple val(sampleId), val(meta), path("${sampleId}.dual_tag_artifact_filter.cutadapt.json"), path("${sampleId}.dual_tag_artifact_filter.summary.tsv"), emit: qc
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def runtimeExports = runtimeShellExports(meta)

    """
    ${runtimeExports}

    "\$PYTHON3_BIN" "${projectDir}/bin/run_dual_tag_artifact_filter.py" \\
      --mode "${mode}" \\
      --sample-id "${sampleId}" \\
      --tagmentation "${meta.dna_tagmentation}" \\
      --r1 "${trimmedR1}" \\
      --r2 "${trimmedR2}" \\
      --signature-fasta "${signatureFasta}" \\
      --cutadapt-bin "\$CUTADAPT_BIN" \\
      --cpus "${task.cpus}" \\
      --output-r1 "${sampleId}.dna_dual_tag_clean.R1.fastq" \\
      --output-r2 "${sampleId}.dna_dual_tag_clean.R2.fastq" \\
      --output-json "${sampleId}.dual_tag_artifact_filter.cutadapt.json" \\
      --output-summary "${sampleId}.dual_tag_artifact_filter.summary.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      cutadapt: "\$("\$CUTADAPT_BIN" --version)"
      component: "local"
    END_VERSIONS
    """
}
