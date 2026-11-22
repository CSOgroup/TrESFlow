include { runtimeShellExports; runtimeOutdir; shellQuoteForFastqCleanup } from '../runtime_support/main'

process BULK_DNA_TARGET {
    tag "${splitName}"
    label 'process_single'
    publishDir { "${runtimeOutdir()}/bulk/dna/${meta.id}/${meta.group}/${meta.mark}" }, mode: 'copy'

    input:
    tuple val(splitName), val(meta), path(bam), path(metrics), path(retention), path(reference, stageAs: 'reference/*')

    output:
    path("${splitName}.cpm.bw")
    path("${splitName}.peaks.bed")
    path("${splitName}.dna.qc.json"), emit: qc

    script:
    def referenceArg = reference ? "--reference ${shellQuoteForFastqCleanup(reference)}" : ''
    """
    ${runtimeShellExports(meta)}
    "\$SAMTOOLS_BIN" sort -N -@ "${task.cpus}" -o by_name.bam "${bam}"
    "\$PYTHON3_BIN" "${projectDir}/bin/bulk_dna.py" "${splitName}" by_name.bam \
        "${meta.canonical_chrom_sizes}" "${retention}" "${metrics}" \
        ${meta.dna_effective_genome_size} ${task.cpus} ${referenceArg} \
        --coverage "\$BAMCOVERAGE_BIN" --epic2 "\$RUNTIME_BIN_DIR/epic2"
    rm by_name.bam
    """
}
