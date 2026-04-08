/*
 * Module: RNA_STARSOLO_ALIGN
 * Runtime command:
 *   bash scripts/core_runtime/RNA_STARSOLO_ALIGN.sh \
 *     <split_name> <tagged.usam> <ref_base_dir> <outdir> <threads> <species>
 *
 * Inputs:
 *   - grouped RNA unmapped SAM from FQ_TO_SAM
 *   - shared RNA reference base dir
 *   - species selector: human | mouse
 * Outputs:
 *   - STARsolo GeneFull directory
 *   - STAR coordinate-sorted aligned BAM used only by the next RNA stage
 */

process RNA_STARSOLO_ALIGN {
    tag "${splitName}"
    label 'codon_wrapper'

    publishDir "${params.outdir ?: "${projectDir}/results"}/align", mode: 'copy', overwrite: true, pattern: "${splitName}.Solo.outGeneFull"

    input:
    tuple val(splitName), val(meta), path(usam), val(refBaseDir), val(species)

    output:
    tuple val(splitName), val(meta), path("${splitName}.Solo.outGeneFull"), emit: solo_dir
    tuple val(splitName), val(meta), path("${splitName}.Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = params.core_scripts_dir ?: "${projectDir}/scripts/core_runtime"

    if( mode == 'mock' ) {
        """
        mkdir -p "${splitName}.Solo.outGeneFull/filtered"

        cat > "${splitName}.Solo.outGeneFull/filtered/barcodes.tsv" <<'EOF'
mock_barcode
EOF

        cat > "${splitName}.Solo.outGeneFull/filtered/features.tsv" <<'EOF'
mock_feature\tmock_feature\tGene Expression
EOF

        cat > "${splitName}.Solo.outGeneFull/filtered/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
1 1 1
1 1 1
EOF

        printf 'mock aligned bam for %s\n' "${splitName}" > "${splitName}.Aligned.sortedByCoord.out.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        if [[ ! -x "\$STAR_BIN" ]]; then
          echo "Missing configured RNA runtime executable: \$STAR_BIN" >&2
          exit 1
        fi

        bash "${coreScriptsDir}/RNA_STARSOLO_ALIGN.sh" \\
          "${splitName}" \\
          "${usam}" \\
          "${refBaseDir}" \\
          "." \\
          "${task.cpus}" \\
          "${species}"
 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
