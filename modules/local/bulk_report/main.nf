include { runtimeShellExports; runtimeOutdir } from '../runtime_support/main'

process BULK_REPORT {
    label 'process_single'
    publishDir { "${runtimeOutdir()}/alignments" }, mode: 'copy'

    input:
    tuple val(meta), val(spec), path(files, stageAs: 'inputs/*')

    output:
    path('bulk_report.html')
    path('gene_by_group.tsv'), optional: true

    script:
    def encoded = groovy.json.JsonOutput.toJson(spec).bytes.encodeBase64().toString()
    """
    ${runtimeShellExports(meta)}
    "\$PYTHON3_BIN" "${projectDir}/bin/bulk_report.py" '${encoded}' inputs
    """
}
