/*
 * Module: TRES_REPORT_HTML
 *
 * Builds a compact TrESFlow-specific HTML report from existing pipeline QC
 * artifacts. This is intentionally separate from MultiQC because TrESFlow has
 * assay-specific RNA/DNA barcode and mapping semantics that should be explained
 * in one stable end-of-run page.
 */

process TRES_REPORT_HTML {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir ?: "${projectDir}/results"}/tres_report", mode: 'copy', overwrite: true, pattern: "tres_report*"

    input:
    tuple val(meta), path(reportInputs, stageAs: "inputs/?/*")

    output:
    tuple val(meta), path("tres_report.html"), emit: html
    tuple val(meta), path("tres_report_metrics.json"), emit: metrics_json
    path("versions.yml"), emit: versions

    script:
    """
    python3 "${projectDir}/bin/render_tres_report.py" \\
      --input-dir inputs \\
      --output-html tres_report.html \\
      --output-json tres_report_metrics.json \\
      --library-name "${meta.library_name ?: 'unknown library'}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
