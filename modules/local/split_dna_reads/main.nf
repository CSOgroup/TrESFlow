/*
 * Module: SPLIT_DNA_READS
 * Upstream reference:
 *   codon run -plugin seq -release Split_ReadsV2.codon \
 *     <Sample> <OutFolder> <LibName> dna <mo_map.tsv> <trimmed_R1.fq.gz> <trimmed_R2.fq.gz> <sb_group_map.tsv>
 *
 * Inputs:
 *   - sample metadata
 *   - trim_galore DNA FASTQs from the CB-tagged reads
 *   - DNA modality map TSV keyed by sample, group, mark, and modality barcode
 *   - shared sample-barcode group map TSV keyed by sample and group
 * Outputs:
 *   - per-group per-mark DNA FASTQ pairs named as final pigz-compressed split outputs
 *   - per-group per-mark SAM RG header TSVs named as upstream Split_ReadsV2 outputs
 */

import RuntimeSupport

process SPLIT_DNA_READS {
    tag "${sampleId}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/dna_split_fastqs", mode: 'copy', overwrite: true, pattern: "${sampleId}_*_R[12].fastq.gz"

    input:
    tuple val(sampleId), val(meta), path(trimmedR1), path(trimmedR2), path(moMap), path(sbGroupMap)

    output:
    tuple val(sampleId), val(meta), path("${sampleId}_*_R1.fastq.gz"), path("${sampleId}_*_R2.fastq.gz"), emit: split_fastqs
    tuple val(sampleId), val(meta), path("SAM_RG_Header_${sampleId}_*.tsv"), emit: rg_headers
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')
    def runtimeExports = RuntimeSupport.shellExports(meta)

    """
    ${runtimeExports}

    "\$PYTHON3_BIN" "${projectDir}/bin/run_split_reads_dna.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/Split_ReadsV2.codon" \\
      --r1 "${trimmedR1}" \\
      --r2 "${trimmedR2}" \\
      --mo-map "${moMap}" \\
      --sb-group-map "${sbGroupMap}" \\
      --sample "${sampleId}" \\
      --library-name "${meta.library_name}" \\
      --output-dir "." \\
      --pigz-threads "${task.cpus}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
