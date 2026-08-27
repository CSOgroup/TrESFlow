/*
 * Metrics-only exact intersection of retained FASTQ pairs with the barcode
 * decisions already recorded by Tag_Lig3. No reads are written or rerouted.
 */

include { runtimeOutdir } from '../runtime_support/main'

process BARCODE_GATE_METRICS {
    tag "${sampleId}:${modality}"
    label 'process_single'

    conda "${moduleDir}/../python_helpers/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://docker.io/library/python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'
        : 'docker.io/library/python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'}"

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*barcode_gates.tsv"
    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*barcode_composition.tsv"

    input:
    tuple val(sampleId), val(meta), val(modality), path(r1), path(r2), path(tagRecords), path(contractMaps)
    path helperScripts, stageAs: "helper/bin/*"

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.${modality}_barcode_gates.tsv"), path("${sampleId}.${modality}_barcode_composition.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def sbGroupMap = contractMaps[0]
    def moArgument = modality == 'dna' ? "--mo-map \"${contractMaps[1]}\"" : ''

    """
    python3 "helper/bin/write_barcode_gate_metrics.py" \\
      --sample "${sampleId}" \\
      --modality "${modality}" \\
      --r1 "${r1}" \\
      --r2 "${r2}" \\
      --tag-records "${tagRecords}" \\
      --sb-group-map "${sbGroupMap}" \\
      ${moArgument} \\
      --output-gates "${sampleId}.${modality}_barcode_gates.tsv" \\
      --output-composition "${sampleId}.${modality}_barcode_composition.tsv"

    printf '%s\n' \\
      '"${task.process}":' \\
      '  python: "'"\$(python3 -c 'import platform; print(platform.python_version())')"'"' \\
      '  component: "local"' \\
      > versions.yml
    """
}
