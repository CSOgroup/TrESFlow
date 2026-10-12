/* Preprocessing-only RNA conversion and gzip publication in one task. */
include { runtimeShellExports; runtimeOutdir; intermediateFastqCleanupCommand; runtimeWorkDir } from '../runtime_support/main'

process EXPORT_RNA_FASTQS {
    tag "${sampleId}"
    label 'process_compression'

    publishDir { "${runtimeOutdir()}/rna_split_fastqs" }, mode: 'copy', overwrite: true, pattern: '*_R[12].fastq.gz'
    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: '*.rna_export_*.tsv'

    input:
    tuple val(sampleId), val(meta), path(splitR1s), path(splitR2s)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fastq.gz"), path("${sampleId}_*_R2.fastq.gz"), emit: fastqs
    tuple val(sampleId), val(meta), path("${sampleId}.rna_export_retention.tsv"), emit: retention_metrics
    path("${sampleId}.rna_export_pools.tsv"), emit: pools
    path('versions.yml'), emit: versions

    script:
    def r1s = splitR1s instanceof List ? splitR1s : [splitR1s]
    def r2s = splitR2s instanceof List ? splitR2s : [splitR2s]
    def r2BySplit = r2s.collectEntries { path -> [(path.getName().replaceFirst(/_R2\.fastq$/, '')): path] }
    if( r1s.size() != r2s.size() || r2BySplit.size() != r2s.size() ) error 'RNA export requires one R1/R2 per split'
    def manifest = r1s.collect { path ->
        def splitName = path.getName().replaceFirst(/_R1\.fastq$/, '')
        if( !r2BySplit.containsKey(splitName) ) error "Missing RNA export mate for ${splitName}"
        [split_name: splitName, group: meta.split_targets[splitName].group, r1: path.toString(), r2: r2BySplit[splitName].toString()]
    }
    def manifestBase64 = groovy.json.JsonOutput.toJson(manifest).bytes.encodeBase64().toString()
    // The current layout is SB[0:4], UM[4:14]. Neither tagging nor paired
    // trimming clips the 5' prefix; this export is its sole removal point.
    def prefixLength = Math.max((meta.sample_bc_start as int) + (meta.sample_bc_len as int), (meta.umi_bc_start as int) + (meta.umi_bc_len as int))
    def runtimeExports = runtimeShellExports(meta)
    def stagedCleanup = intermediateFastqCleanupCommand(params.cleanup_work, 'staged-input', runtimeWorkDir(), "${projectDir}/bin/cleanup_intermediate_fastqs.py", r1s + r2s)
    """
    ${runtimeExports}
    printf '%s' '${manifestBase64}' | base64 --decode > export_manifest.json
    "\$PYTHON3_BIN" "${projectDir}/bin/export_rna_fastqs.py" \\
      --manifest export_manifest.json --sample "${sampleId}" \\
      --prefix-length ${prefixLength} --umi-length ${meta.umi_bc_len} \\
      --minimum-length 20 --pigz "\$PIGZ_BIN" --threads ${task.cpus}
    ${stagedCleanup}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
