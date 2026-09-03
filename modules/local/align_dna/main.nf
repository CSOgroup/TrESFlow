/*
 * Module: ALIGN_DNA
 * Upstream reference:
 *   ./AlignDNA.sh <modality> <sample_name> <R1> <R2> <blacklist.bed> <SAM_RG_Header.tsv> <bwa_prefix> <effective_genome_size> <outdir>
 *
 * Inputs:
 *   - one uncompressed split DNA FASTQ pair from Split_ReadsV2 dna mode
 *   - matching SAM RG header TSV from Split_ReadsV2 dna mode
 *   - bwa-mem2 index prefix inferred from references.dna_ref_dir
 *   - explicit blacklist BED path
 *   - explicit effective genome size
 * Outputs:
 *   - filtered aligned BAM emitted directly by AlignDNA.sh
 *   - BAM index emitted directly by AlignDNA.sh
 *
 * Notes:
 *   - Real execution uses the repo-owned core runtime copy under scripts/core_runtime/.
 *   - AlignDNA.sh reads exported thread settings and keeps proper-pair mapped filtering.
 */

include { runtimeOutdir } from '../runtime_support/main'

process ALIGN_DNA {
    tag "${splitName}"
    label 'dna_processing'

    conda "${moduleDir}/../dna_processing/environment-align.yml"
    container 'community.wave.seqera.io/library/bwa-mem2_samtools@sha256:ce8cbf5cc21c690c8c2994d9bbb409b9313c47f38c5452a5fe7dec3402eff9c8'

    publishDir { "${runtimeOutdir()}/TrES_Stats" }, mode: 'copy', overwrite: true, pattern: "*.dna_alignment_retention.tsv"

    input:
    tuple val(splitName), val(meta), val(sampleGroup), val(modality), path(splitR1), path(splitR2), path(rgHeader), val(bwaReferenceName), path(bwaIndexFiles, stageAs: 'bwa_index/*'), path(blacklistBed), val(effectiveGenomeSize)
    path runtimeScripts, stageAs: 'tresflow/runtime/*'

    output:
    tuple val(splitName), val(meta), path("${splitName}.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}.dna_alignment_retention.tsv"), emit: retention_metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def alignThreads = task.cpus as int
    def viewThreads = alignThreads
    def sortThreads = alignThreads
    def coreScriptsDir = 'tresflow/runtime'
    def bwaReference = "bwa_index/${bwaReferenceName}"

    if( mode == 'mock' ) {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        printf 'mock bam for %s\n' "${splitName}" > "${splitName}.bam"
        printf 'mock bai for %s\n' "${splitName}" > "${splitName}.bam.bai"
        pair_count="\$(awk 'END { print int(NR / 4) }' "${splitR1}")"
        cat > "${splitName}.dna_alignment_retention.tsv" <<EOF
split_id	metric	pairs	unit
${splitName}	bwa_primary_pairs	\${pair_count}	primary_read1_pair_representatives
${splitName}	post_blacklist_primary_pairs	\${pair_count}	primary_read1_pair_representatives
${splitName}	post_blacklist_mapped_primary_pairs	\${pair_count}	primary_read1_pair_representatives
${splitName}	proper_pair_primary_pairs	\${pair_count}	primary_read1_pair_representatives
EOF

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        export ALIGN_DNA_THREADS="${alignThreads}"
        export ALIGN_DNA_VIEW_THREADS="${viewThreads}"
        export ALIGN_DNA_SORT_THREADS="${sortThreads}"
        export ALIGN_DNA_SORT_MEM="1G"

        bash "${coreScriptsDir}/AlignDNA.sh" \\
          "${modality}" \\
          "${sampleGroup}" \\
          "${splitR1}" \\
          "${splitR2}" \\
          "${blacklistBed}" \\
          "${rgHeader}" \\
          "${bwaReference}" \\
          "${effectiveGenomeSize}" \\
          "."

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
