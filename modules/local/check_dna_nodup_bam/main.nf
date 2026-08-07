/*
 * Guard nf-core/deeptools/bamcoverage from empty DNA NoDup BAMs.
 *
 * deepTools can be noisy or unhelpful on BAMs with zero mapped reads. The old
 * local BAM_COVERAGE_DNA module skipped those BAMs and wrote a warning artifact;
 * this module preserves that behavior before routing non-empty BAMs to nf-core.
 */

include { runtimeShellExports; runtimeOutdir } from '../runtime_support/main'

process CHECK_DNA_NODUP_BAM {
    tag "${splitName}"
    label 'process_single'

    publishDir { "${runtimeOutdir()}/pipeline_info/warnings" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*.zero_mapped_nodup_bam.tsv"

    input:
    tuple val(splitName), val(meta), path(noDupBam, stageAs: "input_NoDup.bam"), path(noDupBai, stageAs: "input_NoDup.bam.bai"), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("*_NoDup.bam"), path("*_NoDup.bam.bai"), val(effectiveGenomeSize), optional: true, emit: ready
    tuple val(splitName), val(meta), path("${splitName}.zero_mapped_nodup_bam.tsv"), optional: true, emit: warnings
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def runtimeExports = runtimeShellExports(meta)
    def sampleId = meta.id as String
    def suffix = splitName.replaceFirst("^${sampleId}_", '')
    def tokens = suffix.tokenize('_')
    def groupName = tokens ? tokens[0] : ''
    def markName = tokens.size() > 1 ? tokens[1..-1].join('_') : ''

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        touch "${splitName}_NoDup.bam" "${splitName}_NoDup.bam.bai"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
    """
    ${runtimeExports}

    if [[ ! -x "\$SAMTOOLS_BIN" ]]; then
      echo "Missing configured DNA runtime executable: \$SAMTOOLS_BIN" >&2
      exit 1
    fi

    mapped_reads="\$("\$SAMTOOLS_BIN" view --threads "${task.cpus}" -c -F 4 "${noDupBam}")"
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
      echo "Skipped nf-core/deeptools/bamcoverage for ${splitName}" >&2
      printf 'sample\tgroup\tmark\tbam\tmapped_reads\tskipped_output\n' > "${splitName}.zero_mapped_nodup_bam.tsv"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \\
        "${sampleId}" \\
        "${groupName}" \\
        "${markName}" \\
        "\${bam_path}" \\
        "\${mapped_reads}" \\
        "${splitName}_NoDup.bw" \\
        >> "${splitName}.zero_mapped_nodup_bam.tsv"
    else
      cp -L "${noDupBam}" "${splitName}_NoDup.bam"
      cp -L "${noDupBai}" "${splitName}_NoDup.bam.bai"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
    }
}
