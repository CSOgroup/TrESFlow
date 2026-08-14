/*
 * Metrics-only exact intersection of retained FASTQ pairs with the barcode
 * decisions already recorded by Tag_Lig3. No reads are written or rerouted.
 */

include { runtimeShellExports; runtimeOutdir } from '../runtime_support/main'

process BARCODE_GATE_METRICS {
    tag "${sampleId}:${modality}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*barcode_gates.tsv"
    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*barcode_composition.tsv"

    input:
    tuple val(sampleId), val(meta), val(modality), path(r1), path(r2), path(tagRecords), path(contractMaps)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}.${modality}_barcode_gates.tsv"), path("${sampleId}.${modality}_barcode_composition.tsv"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def runtimeExports = runtimeShellExports(meta)
    def sbGroupMap = contractMaps[0]
    def moArgument = modality == 'dna' ? "--mo-map \"${contractMaps[1]}\"" : ''

    """
    ${runtimeExports}

    "\$PYTHON3_BIN" "${projectDir}/bin/write_barcode_gate_metrics.py" \\
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
      '  component: "local"' \\
      > versions.yml
    """
}
