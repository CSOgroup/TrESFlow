/*
 * Module: SPLIT_DUPLICATES_DNA
 * Upstream reference:
 *   samtools view --threads <threads> --bam --with-header \
 *     --exclude-flags 0x400 \
 *     --output <sample>_NoDup.bam \
 *     <sample>_MarkedDup.bam
 *
 *   samtools index --threads <threads> --bai \
 *     --output <sample>_NoDup.bam.bai \
 *     <sample>_NoDup.bam
 *
 * Inputs:
 *   - canonical, coordinate-sorted duplicate-marked DNA BAM from
 *     NORMALIZE_DNA_MARKDUPLICATES
 *   - effective genome size propagated to the coverage branch
 * Outputs:
 *   - duplicate-filtered NoDup BAM
 *   - NoDup BAM index
 *   - mapped-read count used to gate bamCoverage
 *   - the historical zero-mapped warning when coverage must be skipped
 *
 * Notes:
 *   - No process consumes a duplicate-only BAM, so duplicates are discarded without encoding a temporary file.
 *   - The NoDup definition remains all records without flag 0x400.
 *   - Canonical filtering is intentionally not repeated here: the sole input
 *     producer, NORMALIZE_DNA_MARKDUPLICATES, already enforces that invariant.
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir } from '../runtime_support/main'

process SPLIT_DUPLICATES_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir { "${runtimeOutdir()}/dna_align" }, mode: 'copy', overwrite: true, pattern: "*_NoDup.bam*"
    publishDir { "${runtimeOutdir()}/pipeline_info/warnings" }, mode: params.publish_dir_mode, overwrite: true, pattern: "*.zero_mapped_nodup_bam.tsv"

    input:
    tuple val(splitName), val(meta), path(markedDupBam), val(effectiveGenomeSize)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam.bai"), emit: bai
    tuple val(splitName), val(meta), val(effectiveGenomeSize), path("${splitName}.nodup_mapped_reads.txt"), emit: mapped_reads
    tuple val(splitName), val(meta), path("${splitName}.zero_mapped_nodup_bam.tsv"), optional: true, emit: warnings
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)
    def sampleId = meta.id as String
    def suffix = splitName.replaceFirst("^${sampleId}_", '')
    def tokens = suffix.tokenize('_')
    def groupName = tokens ? tokens[0] : ''
    def markName = tokens.size() > 1 ? tokens[1..-1].join('_') : ''

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock nodup bam for %s\n' "${splitName}" > "${splitName}_NoDup.bam"
        printf 'mock nodup bai for %s\n' "${splitName}" > "${splitName}_NoDup.bam.bai"
        printf '1\n' > "${splitName}.nodup_mapped_reads.txt"

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

        echo "Using SAMTOOLS_BIN=\$SAMTOOLS_BIN"

        bash "${coreScriptsDir}/SplitDuplicatesDNA.sh" \\
          "${markedDupBam}" \\
          "${splitName}_NoDup.bam" \\
          "${splitName}_NoDup.bam.bai" \\
          "${splitName}.nodup_mapped_reads.txt" \\
          "${splitName}.zero_mapped_nodup_bam.tsv" \\
          "${task.cpus}" \\
          "${sampleId}" \\
          "${groupName}" \\
          "${markName}" \\
          "${splitName}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
