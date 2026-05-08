/*
 * Module: SEQUENCING_EFFICIENCY
 * Inputs:
 *   - RNA/DNA tag-record tables emitted by the barcode-tagging steps
 *   - RNA filtered-cell BAMs
 *   - DNA duplicate-marked and NoDup BAMs
 *   - derived SB group and DNA modality maps
 * Outputs:
 *   - sequencing-efficiency count tables, Sankey plots, UpSet plots, combined summaries, and warnings
 */

import RuntimeSupport

process SEQUENCING_EFFICIENCY {
    tag "sequencing_efficiency"

    publishDir "${params.outdir ?: "${projectDir}/results"}/TrES_Stats", mode: 'copy', overwrite: true

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
    path("*.*_sequencing_efficiency*"), emit: reports
    path("sequencing_efficiency.warnings.tsv"), emit: warnings
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
    mkdir -p "\$MPLCONFIGDIR"

    "\$PYTHON3_BIN" "${projectDir}/bin/run_sequencing_efficiency.py" \\
      --outdir . \\
      --rna-tag-records ${rnaTagArgs} \\
      --rna-filtered-bams ${rnaBamArgs} \\
      --dna-tag-records ${dnaTagArgs} \\
      --dna-markeddup-bams ${dnaMarkedDupArgs} \\
      --dna-nodup-bams ${dnaNoDupArgs} \\
      --sb-group-maps ${sbGroupMapArgs} \\
      --dna-mo-maps ${dnaMoMapArgs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      component: "local"
    END_VERSIONS
    """
}
