/*
 * Module: TRES_REPORT_HTML
 *
 * Builds the self-contained TrESFlow QC report and its normalized tables from
 * explicitly staged upstream metric artifacts. The process uses
 * the same repository-owned Python package as the standalone assessor.
 */

include { runtimeOutdir; runtimeShellExports } from '../runtime_support/main'

process TRES_REPORT_HTML {
    tag "${meta.id}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "{tres_report.html,read_retention.tsv,qc_metrics.tsv,barcode_composition.tsv,library_complexity.tsv}"

    input:
    tuple val(meta), path(reportInputs, stageAs: "inputs/?/*")

    output:
    tuple val(meta), path("tres_report.html"), emit: html
    tuple val(meta), path("read_retention.tsv"), emit: read_retention
    tuple val(meta), path("qc_metrics.tsv"), emit: qc_metrics
    tuple val(meta), path("barcode_composition.tsv"), emit: barcode_composition
    tuple val(meta), path("library_complexity.tsv"), emit: library_complexity

    script:
    def reportTitle = meta.report_title?.toString()?.trim()
    def pipelineVersion = meta.pipeline_version?.toString()?.trim()
    if( !reportTitle ) {
        error 'TRES_REPORT_HTML requires the resolved root output-directory basename as meta.report_title'
    }
    if( !pipelineVersion ) {
        error 'TRES_REPORT_HTML requires the launch-resolved repository release version as meta.pipeline_version'
    }
    def metadataJson = groovy.json.JsonOutput.toJson(meta.findAll { key, _value -> !key.toString().startsWith('runtime_') })
    def metadataBase64 = metadataJson.bytes.encodeBase64().toString()
    def runtimeExports = runtimeShellExports(meta)
    """
    ${runtimeExports}
    "\$PYTHON3_BIN" "${projectDir}/bin/render_tres_report.py" \
      --input-dir inputs \
      --output-dir . \
      --output-html tres_report.html \
      --no-json \
      --no-standalone-figures \
      --title "${reportTitle}" \
      --pipeline-version "${pipelineVersion}" \
      --run-metadata-base64 "${metadataBase64}"
    """
}
