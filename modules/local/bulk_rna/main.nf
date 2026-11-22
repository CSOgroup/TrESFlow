include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir; shellQuoteForFastqCleanup } from '../runtime_support/main'

process PREPARE_BULK_RNA {
    tag "${splitName}"
    label 'process_single'

    input:
    tuple val(splitName), val(meta), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*')

    output:
    tuple val(splitName), val(meta), path('input.sam'), path('preparation.json'), emit: prepared

    script:
    def spec = [
        modality: 'rna', sample: meta.id, group: meta.group, fastq_only: true,
        ligation_length: 3 * meta.cell_bc_len, umi_start: meta.umi_bc_start, umi_length: meta.umi_bc_len,
        sb_start: meta.sample_bc_start, sb_length: meta.sample_bc_len, sb_hd: meta.sample_hd,
        sb_rc: meta.sample_reverse_complement == 'rev',
    ] + (meta.bulk_validation ?: [:])
    def encoded = groovy.json.JsonOutput.toJson(spec).bytes.encodeBase64().toString()
    """
    ${runtimeShellExports(meta)}
    "\$PYTHON3_BIN" "${projectDir}/bin/prepare_bulk_rna.py" '${encoded}' \
        ${shellQuoteForFastqCleanup(r1)} ${shellQuoteForFastqCleanup(r2)} input.sam preparation.json
    """
}

process BULK_RNA_MOLECULES {
    tag "${splitName}"
    label 'process_single'
    publishDir { "${runtimeOutdir()}/bulk/rna/${meta.id}/${meta.group}" }, mode: 'copy', pattern: "*.{bam,bai,bw,json}"

    input:
    tuple val(splitName), val(meta), path(aligned), path(solo), path(preparation)

    output:
    path("${splitName}.dedup.bam")
    path("${splitName}.dedup.bam.bai")
    path("${splitName}.rpm.bw")
    path("${splitName}.counts.tsv"), emit: counts
    path("${splitName}.rna.qc.json"), emit: qc

    script:
    """
    ${runtimeShellExports(meta)}
    "\$SAMTOOLS_BIN" sort -N -@ ${task.cpus} -o by_name.bam "${aligned}"
    "\$PYTHON3_BIN" "${projectDir}/bin/bulk_rna.py" "${splitName}" "${meta.rna_star_index_dir}" \
        "${preparation}" by_name.bam "${solo}"
    "\$SAMTOOLS_BIN" sort -@ ${task.cpus} -o "${splitName}.dedup.bam" representatives.bam
    "\$SAMTOOLS_BIN" index "${splitName}.dedup.bam"
    bash "${runtimeCoreScriptsDir()}/RNA_COVERAGE.sh" "${splitName}" "${splitName}.dedup.bam" \
        "${meta.rna_star_index_dir}" "${meta.rna_star_index_dir}/chrNameLength.txt" . "${task.cpus}" bulk
    rm molecules.sqlite representatives.bam by_name.bam
    """
}
