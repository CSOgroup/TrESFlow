/*
 * Module: RNA_QC
 * Runtime command:
 *   python3 bin/run_core_qc.py rna \
 *     --sb-group-map <sb_group_map.tsv> \
 *     --sample-counts <sample_barcode.counts.tsv>... \
 *     --cell-counts <cell.counts.tsv>... \
 *     --solo-dir <sample_group.Solo.outGeneFull>... \
 *     --outdir .
 *
 * Inputs:
 *   - derived sample-barcode group map
 *   - published RNA sample-barcode counts tables
 *   - published RNA cell-barcode counts tables
 *   - published STARsolo GeneFull directories
 * Outputs:
 *   - RNA QC stage-count tables
 *   - per-sample and per-group RNA QC plots under rna/
 */

process RNA_QC {
    tag "rna_qc"

    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    path(sbGroupMap)
    path(sampleCounts)
    path(cellCounts)
    path(soloDirs)

    output:
    path("rna_sample_stage_counts.tsv"), emit: sample_table
    path("rna_group_stage_counts.tsv"), emit: group_table
    path("rna"), emit: plots

    script:
    def sampleCountArgs = sampleCounts.collect { "--sample-counts '${it}'" }.join(" \\\n      ")
    def cellCountArgs = cellCounts.collect { "--cell-counts '${it}'" }.join(" \\\n      ")
    def soloDirArgs = soloDirs.collect { "--solo-dir '${it}'" }.join(" \\\n      ")

    """
    export MPLCONFIGDIR="\$(mktemp -d /tmp/mplconfig-rna-qc.XXXXXX)"
    trap 'rm -rf "\${MPLCONFIGDIR}"' EXIT

    "\$PYTHON3_BIN" "${projectDir}/bin/run_core_qc.py" rna \\
      --sb-group-map "${sbGroupMap}" \\
      ${sampleCountArgs} \\
      ${cellCountArgs} \\
      ${soloDirArgs} \\
      --outdir "."
    """
}
