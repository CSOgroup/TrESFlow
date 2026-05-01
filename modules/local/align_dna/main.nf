/*
 * Module: ALIGN_DNA
 * Upstream reference:
 *   ./AlignDNA.sh <modality> <sample_name> <R1> <R2> <blacklist.bed> <SAM_RG_Header.tsv> <bwa_prefix> <effective_genome_size> <outdir>
 *
 * Inputs:
 *   - one split DNA FASTQ pair from Split_ReadsV2 dna mode
 *   - matching SAM RG header TSV from Split_ReadsV2 dna mode
 *   - bwa-mem2 index prefix inferred from references.dna_ref_dir
 *   - explicit blacklist BED path
 *   - explicit effective genome size
 * Outputs:
 *   - filtered aligned BAM emitted directly by AlignDNA.sh
 *   - BAM index emitted directly by AlignDNA.sh
 *   - properly paired mapped reads per barcode TSV emitted directly by AlignDNA.sh
 *
 * Notes:
 *   - Real execution uses the repo-owned core runtime copy under scripts/core_runtime/.
 *   - AlignDNA.sh reads exported thread settings and keeps the upstream good-barcode threshold.
 */

import RuntimeSupport

process ALIGN_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/dna_align", mode: 'copy', overwrite: true, pattern: "*_ProperPairedMapped_reads_per_barcode.tsv"

    input:
    tuple val(splitName), val(meta), val(sampleGroup), val(modality), path(splitR1), path(splitR2), path(rgHeader), val(bwaReference), val(blacklistBed), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("${splitName}.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}_ProperPairedMapped_reads_per_barcode.tsv"), emit: barcode_counts
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def alignThreads = task.cpus as int
    def viewThreads = alignThreads
    def sortThreads = alignThreads
    def coreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')
    def runtimeExports = RuntimeSupport.shellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock bam for %s\n' "${splitName}" > "${splitName}.bam"
        printf 'mock bai for %s\n' "${splitName}" > "${splitName}.bam.bai"
        cat > "${splitName}_ProperPairedMapped_reads_per_barcode.tsv" <<'EOF'
1	mock_barcode
EOF

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        ${runtimeExports}

        for required_bin in "\$BWA_MEM2_BIN" "\$SAMTOOLS_BIN"; do
          if [[ ! -x "\${required_bin}" ]]; then
            echo "Missing configured DNA runtime executable: \${required_bin}" >&2
            exit 1
          fi
        done

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
