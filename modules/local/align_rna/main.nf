/*
 * Module: ALIGN_RNA
 * Upstream reference:
 *   ./AlignRNA.sh <sample_name> <tagged.usam> <ref_base_dir> <outdir> <threads> <species>
 *
 * Inputs:
 *   - grouped RNA unmapped SAM from FQ_TO_SAM
 *   - shared RNA STAR reference base dir
 *   - species selector: human | mouse
 * Outputs:
 *   - STARsolo GeneFull directory for the grouped RNA sample
 *   - filtered-cells BAM emitted by AlignRNA.sh
 *   - stranded and unstranded bigWig signal tracks when produced by STAR inputAlignmentsFromBAM
 *
 * Notes:
 *   - AlignRNA.sh derives the actual reference paths internally:
 *       human -> <ref_base_dir>/GRCh38_TrES/star and <ref_base_dir>/hg38.chrom.sizes
 *       mouse -> <ref_base_dir>/GRCm39_TrES/star and <ref_base_dir>/mm39.chrom.sizes
 *   - The script auto-detects CB length from CB:Z: or CR:Z: in the input unmapped SAM.
 *   - Real execution uses the repo-owned core runtime copy under scripts/core_runtime/.
 */

process ALIGN_RNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/align", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(usam), val(refBaseDir), val(species)

    output:
    tuple val(splitName), val(meta), path("${splitName}.Solo.outGeneFull"), emit: solo_dir
    tuple val(splitName), val(meta), path("${splitName}.filtered_cells.bam"), emit: filtered_bam
    tuple val(splitName), val(meta), path("${splitName}.stranded_*.bw"), optional: true, emit: stranded_bw
    tuple val(splitName), val(meta), path("${splitName}.unstranded_*.bw"), optional: true, emit: unstranded_bw

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        mkdir -p "${splitName}.Solo.outGeneFull/filtered"

        cat > "${splitName}.Solo.outGeneFull/filtered/barcodes.tsv" <<'EOF'
mock_barcode
EOF

        cat > "${splitName}.Solo.outGeneFull/filtered/features.tsv" <<'EOF'
mock_feature\tmock_feature\tGene Expression
EOF

        cat > "${splitName}.Solo.outGeneFull/filtered/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
1 1 1
1 1 1
EOF

        printf 'mock bam for %s\n' "${splitName}" > "${splitName}.filtered_cells.bam"
        printf 'mock stranded bigwig\n' > "${splitName}.stranded_Signal.Unique.str1.out.bw"
        printf 'mock unstranded bigwig\n' > "${splitName}.unstranded_Signal.Unique.str1.out.bw"
        """
    }
    else {
        """
        for required_bin in "${params.runtime_star}" "${params.runtime_samtools}" "${params.runtime_bedgraph_to_bigwig}"; do
          if [[ ! -x "\${required_bin}" ]]; then
            echo "Missing configured RNA runtime executable: \${required_bin}" >&2
            exit 1
          fi
        done

        export STAR_BIN="${params.runtime_star}"
        export SAMTOOLS_BIN="${params.runtime_samtools}"
        export BEDGRAPH_TO_BIGWIG_BIN="${params.runtime_bedgraph_to_bigwig}"

        bash "${params.core_scripts_dir}/AlignRNA.sh" \\
          "${splitName}" \\
          "${usam}" \\
          "${refBaseDir}" \\
          "." \\
          "${task.cpus}" \\
          "${species}"
        """
    }
}
