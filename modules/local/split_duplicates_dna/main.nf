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
 *   - duplicate-marked DNA BAM from GATK MarkDuplicates
 * Outputs:
 *   - duplicate-filtered NoDup BAM
 *   - NoDup BAM index
 *
 * Notes:
 *   - No process consumes a duplicate-only BAM, so duplicates are discarded without encoding a temporary file.
 *   - The NoDup definition remains all records without flag 0x400.
 */

include { runtimeShellExports; runtimeOutdir; runtimeCoreScriptsDir } from '../runtime_support/main'

process SPLIT_DUPLICATES_DNA {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir { "${runtimeOutdir()}/dna_align" }, mode: 'copy', overwrite: true, pattern: "*_NoDup.bam*"

    input:
    tuple val(splitName), val(meta), path(markedDupBam), path(canonicalChromosomes)

    output:
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam"), emit: bam
    tuple val(splitName), val(meta), path("${splitName}_NoDup.bam.bai"), emit: bai
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)

    if( mode == 'mock' ) {
        """
        ${runtimeExports}

        printf 'mock nodup bam for %s\n' "${splitName}" > "${splitName}_NoDup.bam"
        printf 'mock nodup bai for %s\n' "${splitName}" > "${splitName}_NoDup.bam.bai"

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

        bash "${coreScriptsDir}/FilterCanonicalBam.sh" \\
          "${markedDupBam}" \\
          "${splitName}_NoDup.bam" \\
          "${canonicalChromosomes}" \\
          "${task.cpus}" \\
          normal \\
          --exclude-flags 0x400

        "\$SAMTOOLS_BIN" index \\
          --threads "${task.cpus}" \\
          --bai \\
          --output "${splitName}_NoDup.bam.bai" \\
          "${splitName}_NoDup.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
