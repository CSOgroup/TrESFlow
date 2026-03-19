/*
 * Module: ALIGN_DNA
 * Upstream reference:
 *   ./AlignDNA.sh <modality> <sample_name> <R1> <R2> <blacklist.bed> <SAM_RG_Header.tsv> <bwa_prefix> <effective_genome_size> <outdir>
 *
 * Inputs:
 *   - one split DNA FASTQ pair from Split_ReadsV2 dna mode
 *   - matching SAM RG header TSV from Split_ReadsV2 dna mode
 *   - explicit bwa-mem2 index prefix
 *   - explicit blacklist BED path
 *   - explicit effective genome size
 * Outputs:
 *   - filtered aligned BAM emitted directly by AlignDNA.sh
 *   - BAM index emitted directly by AlignDNA.sh
 *   - properly paired mapped reads per barcode TSV emitted directly by AlignDNA.sh
 *
 * Notes:
 *   - This wrapper preserves the upstream shell script as the real execution path.
 *   - AlignDNA.sh internally hardcodes threads and the good-barcode threshold.
 */

process ALIGN_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/dna_align", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), val(sampleGroup), val(modality), path(splitR1), path(splitR2), path(rgHeader), val(bwaReference), val(blacklistBed), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("${splitName}.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}_ProperPairedMapped_reads_per_barcode.tsv"), emit: barcode_counts

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock bam for %s\n' "${splitName}" > "${splitName}.bam"
        printf 'mock bai for %s\n' "${splitName}" > "${splitName}.bam.bai"
        cat > "${splitName}_ProperPairedMapped_reads_per_barcode.tsv" <<'EOF'
1	mock_barcode
EOF
        """
    }
    else {
        """
        bash "${params.upstream_dir}/AlignDNA.sh" \\
          "${modality}" \\
          "${sampleGroup}" \\
          "${splitR1}" \\
          "${splitR2}" \\
          "${blacklistBed}" \\
          "${rgHeader}" \\
          "${bwaReference}" \\
          "${effectiveGenomeSize}" \\
          "."
        """
    }
}
