/*
 * When --publish_split_fastqs is enabled, compress final split FASTQ copies
 * for publication without delaying their uncompressed computational branch.
 * The input FASTQs are staged read-only; pigz -c preserves them for FQ_TO_SAM
 * or ALIGN_DNA. The RNA/DNA subworkflows invoke this process only when
 * publication is enabled.
 */

include { runtimeOutdir } from '../runtime_support/main'

process COMPRESS_SPLIT_FASTQS {
    tag "${modality}.${sampleId}"
    label 'process_single'

    conda "${moduleDir}/../fastq_preprocessing/environment-pigz.yml"
    container 'quay.io/biocontainers/cutadapt@sha256:2049f305574854edb189ccad7038fda4801ef16458bcde0239383d42d4a3f83a'

    publishDir { "${runtimeOutdir()}/${modality}_split_fastqs" }, mode: 'copy', overwrite: true, pattern: "*_R[12].fastq.gz"

    input:
    tuple val(sampleId), val(meta), val(modality), path(splitR1s), path(splitR2s)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fastq.gz"), path("${sampleId}_*_R2.fastq.gz"), emit: compressed_fastqs
    path("versions.yml"), emit: versions

    script:
    def fastqs = (splitR1s instanceof List ? splitR1s : [splitR1s]) + (splitR2s instanceof List ? splitR2s : [splitR2s])
    def fastqArgs = fastqs.collect { fastq -> "\"${fastq}\"" }.join(' ')

    """
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    for fastq in ${fastqArgs}; do
      pigz -c -p "${task.cpus}" "\$fastq" > "\$(basename "\$fastq").gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
