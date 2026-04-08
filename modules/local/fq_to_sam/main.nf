/*
 * Module: FQ_TO_SAM
 * Upstream reference:
 *   codon run -plugin seq -release FqToSAM.codon \
 *     <R1.fq[.gz]> <R2.fq[.gz]> <out.sam>
 *
 * Inputs:
 *   - split RNA FASTQ pair from Split_ReadsV2 rna mode
 * Outputs:
 *   - unmapped SAM carrying CR:Z:CB+UM plus preserved non-CB/non-UM tags
 *
 * Notes:
 *   - Upstream MAINLAUNCH.sh comments say the script must accept gzipped FASTQs.
 *     The checked-in FqToSAM.codon does accept `.gz` inputs directly, and this wrapper preserves that contract.
 */

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
    def coreScriptsDir = params.core_scripts_dir ?: "${projectDir}/scripts/core_runtime"

    """
    "\$PYTHON3_BIN" "${projectDir}/bin/run_fq_to_sam.py" \\
      --mode "${mode}" \\
      --script "${coreScriptsDir}/FqToSAM.codon" \\
      --r1 "${splitR1}" \\
      --r2 "${splitR2}" \\
      --output-sam "${splitName}_tagged.usam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
