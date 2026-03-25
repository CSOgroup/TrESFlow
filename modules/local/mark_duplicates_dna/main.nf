/*
 * Module: MARK_DUPLICATES_DNA
 * Upstream reference:
 *   gatk MarkDuplicates \
 *     -I <input.bam> \
 *     -O <sample>_MarkedDup.bam \
 *     -M <sample>.DuplicateMetrics.txt \
 *     --REMOVE_DUPLICATES false \
 *     --BARCODE_TAG CB \
 *     --CREATE_INDEX true \
 *     --MAX_RECORDS_IN_RAM 10000000
 *
 * Inputs:
 *   - aligned DNA BAM from AlignDNA.sh
 * Outputs:
 *   - duplicate-marked BAM
 *   - BAM index
 *   - duplicate metrics TSV
 *
 * Notes:
 *   - This wrapper preserves the launcher's GATK MarkDuplicates invocation.
 *   - It normalizes the created BAM index to <sample>_MarkedDup.bam.bai for a stable Nextflow output contract.
 */

process MARK_DUPLICATES_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/dna_dedup", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(alignedBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}.DuplicateMetrics.txt"), emit: metrics

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock marked-dup bam for %s\n' "${splitName}" > "${splitName}_MarkedDup.bam"
        printf 'mock marked-dup bai for %s\n' "${splitName}" > "${splitName}_MarkedDup.bam.bai"
        cat > "${splitName}.DuplicateMetrics.txt" <<'EOF'
## mock duplicate metrics
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
mock\t0\t0\t0\t0\t0\t0\t0\t0.0\t0
EOF
        """
    }
    else {
        """
        if [[ ! -x "\$GATK_BIN" ]]; then
          echo "Missing configured GATK executable at \$GATK_BIN" >&2
          exit 1
        fi

        echo "Using GATK_BIN=\$GATK_BIN"

        "\$GATK_BIN" MarkDuplicates \\
          -I "${alignedBam}" \\
          -O "${splitName}_MarkedDup.bam" \\
          -M "${splitName}.DuplicateMetrics.txt" \\
          --REMOVE_DUPLICATES false \\
          --BARCODE_TAG CB \\
          --CREATE_INDEX true \\
          --MAX_RECORDS_IN_RAM 10000000

        if [[ -f "${splitName}_MarkedDup.bai" ]]; then
          mv "${splitName}_MarkedDup.bai" "${splitName}_MarkedDup.bam.bai"
        fi

        if [[ ! -f "${splitName}_MarkedDup.bam.bai" ]]; then
          echo "Missing BAM index for ${splitName}_MarkedDup.bam" >&2
          exit 1
        fi
        """
    }
}
