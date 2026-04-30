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

import RuntimeSupport

process BAM_COVERAGE_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/pipeline_info/warnings", mode: 'copy', overwrite: true, pattern: "*.zero_mapped_nodup_bam.tsv"

    input:
    tuple val(splitName), val(meta), path(noDupBam), path(noDupBai), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bw"), optional: true, emit: bw
    tuple val(splitName), val(meta), path("${splitName}.zero_mapped_nodup_bam.tsv"), optional: true, emit: warnings
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def runtimeExports = RuntimeSupport.shellExports(meta)
    def sampleId = meta.id as String
    def suffix = splitName.replaceFirst("^${sampleId}_", '')
    def tokens = suffix.tokenize('_')
    def groupName = tokens ? tokens[0] : ''
    def markName = tokens.size() > 1 ? tokens[1..-1].join('_') : ''

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock bigwig for %s\n' "${splitName}" > "${splitName}_NoDup.bw"

        printf '%s\n' \
          '"${task.process}":' \
          '  component: "local"' \
          > versions.yml
        """
    }
    else {
        """
        ${runtimeExports}

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
          bam_path="\$(readlink -f "${noDupBam}")"
          cat >&2 <<'EOF'
================================================================================
WARNING: ZERO MAPPED READS IN DNA NoDup BAM
================================================================================
EOF
          echo "Sample: ${sampleId}" >&2
          echo "Group: ${groupName}" >&2
          echo "Mark: ${markName}" >&2
          echo "BAM: \${bam_path}" >&2
          echo "Mapped reads: \${mapped_reads}" >&2
          echo "Skipped output: ${splitName}_NoDup.bw" >&2
          echo "A warning artifact will be published to pipeline_info/warnings/." >&2
          printf 'sample\tgroup\tmark\tbam\tmapped_reads\tskipped_output\n' > "${splitName}.zero_mapped_nodup_bam.tsv"
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${sampleId}" \
            "${groupName}" \
            "${markName}" \
            "\${bam_path}" \
            "\${mapped_reads}" \
            "${splitName}_NoDup.bw" \
            >> "${splitName}.zero_mapped_nodup_bam.tsv"
          printf '%s\n' \\
            '"${task.process}":' \\
            '  component: "local"' \\
            > versions.yml
          exit 0
        fi

        tmp_root="\${TMPDIR:?TMPDIR must be set from samplesheet runtime.tmpdir}"
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
