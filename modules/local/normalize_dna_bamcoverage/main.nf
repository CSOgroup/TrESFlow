/*
 * Normalize nf-core/deeptools/bamcoverage bigWig naming back to TrESFlow's
 * historical <split>_NoDup.bw output contract.
 */

include { runtimeOutdir } from '../runtime_support/main'

process NORMALIZE_DNA_BAMCOVERAGE {
    tag "${splitName}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/dna_align" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*_NoDup.bw"

    input:
    tuple val(splitName), val(meta), path(bigwig)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bw"), emit: bw
    path("versions.yml"), emit: versions

    script:
    """
    if [[ "\$(readlink -f "${bigwig}")" != "\$(readlink -f "${splitName}_NoDup.bw" 2>/dev/null || true)" ]]; then
      cp -L "${bigwig}" "${splitName}_NoDup.bw"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
