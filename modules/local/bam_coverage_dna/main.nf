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

    input:
    tuple val(splitName), val(meta), path(noDupBam), path(noDupBai), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bw"), optional: true, emit: bw
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'

    if( mode == 'mock' ) {
        """
        printf 'mock bigwig for %s\n' "${splitName}" > "${splitName}_NoDup.bw"

        printf '%s\n' \
          '"${task.process}":' \
          '  component: "local"' \
          > versions.yml
        """
    }
    else {
        """
        for required_bin in "\$SAMTOOLS_BIN" "\$BAMCOVERAGE_BIN"; do
          if [[ ! -x "\${required_bin}" ]]; then
            echo "Missing configured DNA runtime executable: \${required_bin}" >&2
            exit 1
          fi
        done

        echo "Using SAMTOOLS_BIN=\$SAMTOOLS_BIN"
        echo "Using BAMCOVERAGE_BIN=\$BAMCOVERAGE_BIN"

        mapped_reads="\$("\$SAMTOOLS_BIN" view -c -F 4 "${noDupBam}")"
        if [[ "\${mapped_reads}" -eq 0 ]]; then
          echo "Skipping bamCoverage for ${splitName}: ${noDupBam} has zero mapped reads" >&2
          printf '%s\n' \\
            '"${task.process}":' \\
            '  component: "local"' \\
            > versions.yml
          exit 0
        fi

        tmp_root="\${TMPDIR:-\$PWD/.tmp}"
        mkdir -p "\${tmp_root}"
        export MPLCONFIGDIR="\$(mktemp -d "\${tmp_root}/mplconfig-${splitName}.XXXXXX")"
        trap 'rm -rf "\${MPLCONFIGDIR}"' EXIT

        "\$BAMCOVERAGE_BIN" \\
          -p "${task.cpus}" \\
          -bs 100 \\
          --extendReads \\
          --centerReads \\
          -b "${noDupBam}" \\
          -o "${splitName}_NoDup.bw" \\
          -of bigwig \\
          --effectiveGenomeSize "${effectiveGenomeSize}"

        printf '%s\n' \\
          '"${task.process}":' \\
          '  component: "local"' \\
          > versions.yml
        """
    }
}
