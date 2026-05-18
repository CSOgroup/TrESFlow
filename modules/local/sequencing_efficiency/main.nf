/*
 * Module: SEQUENCING_EFFICIENCY
 * Inputs:
 *   - RNA/DNA tag-record tables emitted by the barcode-tagging steps
 *   - RNA filtered-cell BAMs
 *   - DNA duplicate-marked and NoDup BAMs
 *   - derived SB group and DNA modality maps
 * Outputs:
 *   - sequencing-efficiency UpSet PDF plots
 */

import RuntimeSupport

process SEQUENCING_EFFICIENCY {
    tag "sequencing_efficiency"

    publishDir "${params.outdir ?: "${projectDir}/results"}/TrES_Stats", mode: 'copy', overwrite: true, pattern: "*.upset.pdf"

    input:
    val runtimeMeta
    path rnaTagRecords
    path rnaFilteredBams
    path dnaTagRecords
    path dnaMarkedDupBams
    path dnaNoDupBams
    path sbGroupMaps
    path dnaMoMaps

    output:
    path("*.upset.pdf"), emit: reports
    path("versions.yml"), emit: versions

    script:
    def runtimeExports = RuntimeSupport.shellExports(runtimeMeta)
    def rnaTagArgs = rnaTagRecords ? rnaTagRecords.collect { it.toString() }.join(' ') : ''
    def rnaBamArgs = rnaFilteredBams ? rnaFilteredBams.collect { it.toString() }.join(' ') : ''
    def dnaTagArgs = dnaTagRecords ? dnaTagRecords.collect { it.toString() }.join(' ') : ''
    def dnaMarkedDupArgs = dnaMarkedDupBams ? dnaMarkedDupBams.collect { it.toString() }.join(' ') : ''
    def dnaNoDupArgs = dnaNoDupBams ? dnaNoDupBams.collect { it.toString() }.join(' ') : ''
    def sbGroupMapArgs = sbGroupMaps ? sbGroupMaps.collect { it.toString() }.join(' ') : ''
    def dnaMoMapArgs = dnaMoMaps ? dnaMoMaps.collect { it.toString() }.join(' ') : ''

    """
    ${runtimeExports}
    export MPLCONFIGDIR="\$TMPDIR/matplotlib"
    export SEQUENCING_EFFICIENCY_TMP="\$TMPDIR/sequencing_efficiency_sort"
    mkdir -p "\$MPLCONFIGDIR" "\$SEQUENCING_EFFICIENCY_TMP"

    "\$PYTHON3_BIN" "${projectDir}/bin/run_sequencing_efficiency.py" \\
      --outdir . \\
      --rna-tag-records ${rnaTagArgs} \\
      --rna-filtered-bams ${rnaBamArgs} \\
      --dna-tag-records ${dnaTagArgs} \\
      --dna-markeddup-bams ${dnaMarkedDupArgs} \\
      --dna-nodup-bams ${dnaNoDupArgs} \\
      --sb-group-maps ${sbGroupMapArgs} \\
      --dna-mo-maps ${dnaMoMapArgs} \\
      --min-read-pairs-per-cell "${params.efficiency_min_read_pairs_per_cell}" \\
      --sort-parallel "${task.cpus}" \\
      --sort-buffer "${params.sequencing_efficiency_sort_buffer}" \\
      --tmpdir "\$SEQUENCING_EFFICIENCY_TMP"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
