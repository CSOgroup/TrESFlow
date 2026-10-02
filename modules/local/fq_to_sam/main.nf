/*
 * Module: FQ_TO_SAM
 * Upstream reference:
 *   codon run -plugin seq -release FqToSAM.codon \
 *     <R1.fq[.gz]> <R2.fq[.gz]> <out.sam>
 *
 * Inputs:
 *   - split RNA FASTQ pair from Split_ReadsV2 rna mode
 * Outputs:
 *   - unmapped SAM carrying full CB/XI and separate molecular UM/UR tags
 *
 * Notes:
 *   - The computational branch supplies plain split FASTQs, avoiding a decode before SAM conversion.
 *   - The checked-in FqToSAM.codon remains compatible with legacy `.gz` inputs.
 */

include { runtimeShellExports; runtimeCoreScriptsDir; intermediateFastqCleanupCommand; runtimeWorkDir } from '../runtime_support/main'

process FQ_TO_SAM {
    tag "${splitName}"
    label 'codon_wrapper'

    input:
    tuple val(splitName), val(meta), path(splitR1), path(splitR2)

    output:
    tuple val(splitName), val(meta), path("${splitName}_tagged.usam"), emit: usam
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = runtimeCoreScriptsDir()
    def runtimeExports = runtimeShellExports(meta)
    def stagedFastqCleanup = intermediateFastqCleanupCommand(params.cleanup_work, 'staged-input', runtimeWorkDir(), "${projectDir}/bin/cleanup_intermediate_fastqs.py", [splitR1, splitR2])

    """
    ${runtimeExports}

    "\$PYTHON3_BIN" "${projectDir}/bin/run_fq_to_sam.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/FqToSAM.codon" \\
      --r1 "${splitR1}" \\
      --r2 "${splitR2}" \\
      --output-sam "${splitName}_tagged.usam"

    ${stagedFastqCleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS

    """
}
