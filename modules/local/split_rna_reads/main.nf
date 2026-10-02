/*
 * Module: SPLIT_RNA_READS
 * Upstream reference:
 *   codon run -plugin seq -release Split_ReadsV2.codon \
 *     <Sample> <OutFolder> <LibName> rna - <trimmed_R1.fq> <trimmed_R2.fq> <sb_group_map.tsv>
 *
 * Inputs:
 *   - sample metadata
 *   - uncompressed trim_galore RNA FASTQs from the CB-tagged reads
 *   - shared sample-barcode group map TSV keyed by sample and group
 * Outputs:
 *   - uncompressed per-group RNA FASTQ pairs for downstream computation
 *   - per-group SAM RG header TSVs named as upstream Split_ReadsV2 outputs
 *
 * The versioned SB map carries explicit normalized indices and upstream injection format.
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir; intermediateFastqCleanupCommand; runtimeWorkDir } from '../runtime_support/main'

process SPLIT_RNA_READS {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.rna_read_retention.tsv"

    input:
    tuple val(sampleId), val(meta), path(trimmedR1), path(trimmedR2), path(sbGroupMap)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fastq"), path("${sampleId}_*_R2.fastq"), emit: split_fastqs
    tuple val(sampleId), val(meta), path("SAM_RG_Header_${sampleId}_*.tsv"), emit: rg_headers
    tuple val(sampleId), val(meta), path("${sampleId}.rna_read_retention.tsv"), emit: retention_metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)
    def stagedFastqCleanup = intermediateFastqCleanupCommand(params.cleanup_work, 'staged-input', runtimeWorkDir(), "${projectDir}/bin/cleanup_intermediate_fastqs.py", [trimmedR1, trimmedR2])

    """
    ${runtimeExports}

    "\$PYTHON3_BIN" "${projectDir}/bin/run_split_reads_rna.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/Split_ReadsV2.codon" \\
      --r1 "${trimmedR1}" \\
      --r2 "${trimmedR2}" \\
      --sb-group-map "${sbGroupMap}" \\
      --sample "${sampleId}" \\
      --library-name "${meta.library_name}" \\
      --output-dir "."

    printf '%s\\n' \\
      '"${task.process}":' \\
      '  component: "local"' \\
      > versions.yml

    ${stagedFastqCleanup}
    """
}
