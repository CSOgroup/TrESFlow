include { runtimeOutdir; runtimeShellExports; shellQuoteForFastqCleanup } from '../runtime_support/main'

process WRITE_PREPROCESSED_MANIFEST {
    label 'process_single'

    publishDir { runtimeOutdir() }, mode: 'copy', pattern: 'restart_manifest.json'
    publishDir { "${runtimeOutdir()}/restart" }, mode: 'copy', pattern: '{preprocessing_qc,dna}'

    input:
    val(meta)
    val(splits)
    path(readGroups, stageAs: 'read_groups/?/*')
    tuple path(metrics, stageAs: 'metrics/?/*'), val(reportIndices), val(multiqcIndices)

    output:
    path 'restart_manifest.json'
    path 'preprocessing_qc'
    path 'dna', optional: true

    script:
    def reportArgs = reportIndices.collect { shellQuoteForFastqCleanup(metrics[it]) }.join(' ')
    def multiqcArgs = multiqcIndices.collect { shellQuoteForFastqCleanup(metrics[it]) }.join(' ')
    def rgArgs = [readGroups].flatten().collect { shellQuoteForFastqCleanup(it) }.join(' ')
    def spec = [
        version: 1, splits: splits, pipeline_version: meta.pipeline_version,
        filter_dual_tag_artifacts: params.filter_dual_tag_artifacts,
    ]
    def options = groovy.json.JsonOutput.toJson(spec).bytes.encodeBase64().toString()
    """
    ${runtimeShellExports(meta)}
    "\$PYTHON3_BIN" "${projectDir}/bin/preprocessed_fastqs.py" write "${options}" \
        --report ${reportArgs} --multiqc ${multiqcArgs} --read-groups ${rgArgs}
    """
}
