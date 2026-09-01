/*
 * Module: RNA_STARSOLO_ALIGN
 * Runtime command:
 *   bash scripts/core_runtime/RNA_STARSOLO_ALIGN.sh \
 *     <split_name> <tagged.usam> <star_index_dir> <outdir> <threads>
 *
 * Inputs:
 *   - grouped RNA unmapped SAM from FQ_TO_SAM
 *   - exact RNA STAR index directory resolved from references.rna_ref_dir
 * Outputs:
 *   - STARsolo GeneFull directory
 *   - STAR coordinate-sorted aligned BAM used only by the next RNA stage
 */

include { runtimeOutdir } from '../runtime_support/main'

process RNA_STARSOLO_ALIGN {
    tag "${splitName}"
    label 'rna_alignment'

    conda "${moduleDir}/../rna_alignment/environment-star.yml"
    container 'quay.io/biocontainers/star@sha256:d5139d9e29bc8a871a020ca5b21bc0599fac428bd061d4b13d863f4cb4400d31'

    publishDir { "${runtimeOutdir()}/rna_align" }, mode: 'copy', overwrite: true, pattern: "*.Solo.outGeneFull"
    publishDir { "${runtimeOutdir()}/rna_align" }, mode: 'copy', overwrite: true, pattern: "*.Log.final.out"

    input:
    tuple val(splitName), val(meta), path(usam), path(starIndexDir)
    path runtimeScripts, stageAs: 'tresflow/runtime/*'

    output:
    tuple val(splitName), val(meta), path("${splitName}.Solo.outGeneFull"), emit: solo_dir
    tuple val(splitName), val(meta), path("${splitName}.Solo.outGeneFull/Summary.csv"), emit: solo_summary
    tuple val(splitName), val(meta), path("report/${splitName}.Solo.summary.csv"), emit: report_solo_summary
    tuple val(splitName), val(meta), path("${splitName}.Aligned.sortedByCoord.out.bam"), emit: aligned_bam
    tuple val(splitName), val(meta), path("${splitName}.Log.final.out"), emit: star_log
    path("versions.yml"), emit: versions

    script:
    def mode = task.ext.mock ? 'mock' : 'real'
    def coreScriptsDir = 'tresflow/runtime'
    if( mode == 'mock' ) {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        input_pairs="\$(awk '\$1 !~ /^@/ { n++ } END { print int(n / 2) }' "${usam}")"
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

        printf 'mock aligned bam for %s %s pairs\n' "${splitName}" "\${input_pairs}" > "${splitName}.Aligned.sortedByCoord.out.bam"
        cat > "${splitName}.Log.final.out" <<EOF
                             Started job on |	mock
                             Started mapping on |	mock
                                    Finished on |	mock
       Mapping speed, Million of reads per hour |	1.00
                          Number of input reads |	\${input_pairs}
                      Average input read length |	100
                                    UNIQUE READS:
                   Uniquely mapped reads number |	\${input_pairs}
                        Uniquely mapped reads % |	100.00%
                          Average mapped length |	100.00
                       Number of splices: Total |	0
            Number of splices: Annotated (sjdb) |	0
                       Number of splices: GT/AG |	0
                       Number of splices: GC/AG |	0
                       Number of splices: AT/AC |	0
               Number of splices: Non-canonical |	0
                      Mismatch rate per base, % |	0.10%
                         Deletion rate per base |	0.00%
                        Deletion average length |	0.00
                        Insertion rate per base |	0.00%
                       Insertion average length |	0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	0
             % of reads mapped to multiple loci |	0.00%
        Number of reads mapped to too many loci |	0
             % of reads mapped to too many loci |	0.00%
                                  UNMAPPED READS:
  Number of reads unmapped: too many mismatches |	0
       % of reads unmapped: too many mismatches |	0.00%
            Number of reads unmapped: too short |	0
                 % of reads unmapped: too short |	0.00%
                Number of reads unmapped: other |	0
                     % of reads unmapped: other |	0.00%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%
EOF
        cat > "${splitName}.Solo.outGeneFull/Summary.csv" <<EOF
Number of Reads,\${input_pairs}
Sequencing Saturation,0.25
Reads Mapped to Genome: Unique+Multiple,1.00
Reads Mapped to Genome: Unique,1.00
Reads Mapped to GeneFull: Unique+Multiple GeneFull,1.00
Reads Mapped to GeneFull: Unique GeneFull,1.00
Estimated Number of Cells,1
UMIs in Cells,100
EOF

        mkdir -p report
        cp "${splitName}.Solo.outGeneFull/Summary.csv" \
          "report/${splitName}.Solo.summary.csv"
        chmod -R a+rX "${splitName}.Solo.outGeneFull"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
    else {
        """
        export TMPDIR="\$PWD/.tmp"
        mkdir -p "\$TMPDIR"

        bash "${coreScriptsDir}/RNA_STARSOLO_ALIGN.sh" \\
          "${splitName}" \\
          "${usam}" \\
          "${starIndexDir}" \\
          "." \\
          "${task.cpus}"

        mkdir -p report
        cp "${splitName}.Solo.outGeneFull/Summary.csv" \
          "report/${splitName}.Solo.summary.csv"

        chmod -R a+rX "${splitName}.Solo.outGeneFull"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
          component: "local"
        END_VERSIONS
        """
    }
}
