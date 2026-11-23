include { runtimeShellExports; runtimeOutdir; runtimeWorkDir; intermediateFastqCleanupCommand; shellQuoteForFastqCleanup } from '../runtime_support/main'

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
    def stagedCleanup = intermediateFastqCleanupCommand(meta.restart_input ? false : params.cleanup_work,
        'staged-input', runtimeWorkDir(), "${projectDir}/bin/cleanup_intermediate_fastqs.py", [r1, r2])
    def encoded = groovy.json.JsonOutput.toJson(spec).bytes.encodeBase64().toString()
    """
    ${runtimeShellExports(meta)}
    "\$PYTHON3_BIN" "${projectDir}/bin/prepare_bulk_rna.py" '${encoded}' \
        ${shellQuoteForFastqCleanup(r1)} ${shellQuoteForFastqCleanup(r2)} input.sam preparation.json
    ${stagedCleanup}
    """
}

process BULK_RNA_MOLECULES {
    tag "${splitName}"
    label 'process_single'
    publishDir { "${runtimeOutdir()}/alignments" }, mode: 'copy', pattern: "*.{bam,bai,bw,json}"

    input:
    tuple val(splitName), val(meta), path(aligned), path(solo), path(starLog), path(preparation)

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
    "\$SAMTOOLS_BIN" collate -u -O "${splitName}.dedup.bam" | \
        "\$PYTHON3_BIN" "${projectDir}/bin/bulk_rna_coverage.py" blocks | \
        LC_ALL=C sort -k1,1n -k2,2n -T . | \
        "\$PYTHON3_BIN" "${projectDir}/bin/bulk_rna_coverage.py" track "${splitName}.dedup.bam" \
        "${splitName}.rna.qc.json" "${splitName}.rpm.bw"
    rm molecules.sqlite representatives.bam by_name.bam
    """
}
