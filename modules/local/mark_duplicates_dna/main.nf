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

import RuntimeSupport

process MARK_DUPLICATES_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/dna_align", mode: 'copy', overwrite: true, pattern: "${splitName}_MarkedDup.bam*"

    input:
    tuple val(splitName), val(meta), path(alignedBam)

    output:
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_MarkedDup.bam.bai"), emit: bai
    tuple val(splitName), val(meta), path("${splitName}.DuplicateMetrics.txt"), emit: metrics
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def runtimeExports = RuntimeSupport.shellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock marked-dup bam for %s\n' "${splitName}" > "${splitName}_MarkedDup.bam"
        printf 'mock marked-dup bai for %s\n' "${splitName}" > "${splitName}_MarkedDup.bam.bai"
        cat > "${splitName}.DuplicateMetrics.txt" <<'EOF'
## mock duplicate metrics
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
mock\t0\t0\t0\t0\t0\t0\t0\t0.0\t0
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
