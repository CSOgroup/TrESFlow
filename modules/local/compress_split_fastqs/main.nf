/*
 * Compress final split FASTQ copies for publication without delaying their
 * uncompressed computational branch. The input FASTQs are staged read-only;
 * pigz -c preserves them for FQ_TO_SAM or ALIGN_DNA.
 */

import RuntimeSupport

process COMPRESS_SPLIT_FASTQS {
    tag "${modality}.${sampleId}"
    label 'process_single'

    publishDir "${params.outdir ?: "${projectDir}/results"}/${modality}_split_fastqs", mode: 'copy', overwrite: true, pattern: "${sampleId}_*_R[12].fastq.gz"

    input:
    tuple val(sampleId), val(meta), val(modality), path(splitR1s), path(splitR2s)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fastq.gz"), path("${sampleId}_*_R2.fastq.gz"), emit: compressed_fastqs
    path("versions.yml"), emit: versions

    script:
    def runtimeExports = RuntimeSupport.shellExports(meta)
    def fastqs = (splitR1s instanceof List ? splitR1s : [splitR1s]) + (splitR2s instanceof List ? splitR2s : [splitR2s])
    def fastqArgs = fastqs.collect { fastq -> "\"${fastq}\"" }.join(' ')

    """
    ${runtimeExports}

    if [[ ! -x "\$PIGZ_BIN" ]]; then
      echo "Missing configured pigz executable for split FASTQ publication: \$PIGZ_BIN" >&2
      exit 1
    fi

    for fastq in ${fastqArgs}; do
      "\$PIGZ_BIN" -c -p "${task.cpus}" "\$fastq" > "\$(basename "\$fastq").gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
