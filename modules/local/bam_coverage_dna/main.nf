/*
 * Module: BAM_COVERAGE_DNA
 * Upstream reference:
 *   bamCoverage -p <threads> -bs 100 --extendReads --centerReads \
 *     -b <sample>_NoDup.bam \
 *     -o <sample>_NoDup.bw \
 *     -of bigwig \
 *     --effectiveGenomeSize <effective_genome_size>
 *
 * Inputs:
 *   - duplicate-filtered DNA NoDup BAM
 *   - matching NoDup BAM index, used to preserve the upstream order
 *   - explicit effective genome size
 * Outputs:
 *   - bigWig coverage track from the NoDup BAM when mapped reads are present
 *
 * Notes:
 *   - This wrapper preserves the launcher's immediate bamCoverage step after NoDup indexing.
 *   - If the NoDup BAM has zero mapped reads, the wrapper skips bamCoverage with a warning.
 *     This is an integration guard for the current real data, where deepTools otherwise floods
 *     stderr and does not finish cleanly on an empty mapped-read set.
 */

process BAM_COVERAGE_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir}/dna_coverage", mode: 'copy', overwrite: true

    input:
    tuple val(splitName), val(meta), path(noDupBam), path(noDupBai), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bw"), optional: true, emit: bw

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock bigwig for %s\n' "${splitName}" > "${splitName}_NoDup.bw"
        """
    }
    else {
        """
        for required_bin in "${params.runtime_samtools}" "${params.runtime_bam_coverage}"; do
          if [[ ! -x "\${required_bin}" ]]; then
            echo "Missing configured DNA runtime executable: \${required_bin}" >&2
            exit 1
          fi
        done

        echo "Using SAMTOOLS_BIN=${params.runtime_samtools}"
        echo "Using BAMCOVERAGE_BIN=${params.runtime_bam_coverage}"

        mapped_reads="\$("${params.runtime_samtools}" view -c -F 4 "${noDupBam}")"
        if [[ "\${mapped_reads}" -eq 0 ]]; then
          echo "Skipping bamCoverage for ${splitName}: ${noDupBam} has zero mapped reads" >&2
          exit 0
        fi

        export MPLCONFIGDIR="\$(mktemp -d /tmp/mplconfig-${splitName}.XXXXXX)"
        trap 'rm -rf "\${MPLCONFIGDIR}"' EXIT

        "${params.runtime_bam_coverage}" \\
          -p "${task.cpus}" \\
          -bs 100 \\
          --extendReads \\
          --centerReads \\
          -b "${noDupBam}" \\
          -o "${splitName}_NoDup.bw" \\
          -of bigwig \\
          --effectiveGenomeSize "${effectiveGenomeSize}"
        """
    }
}
