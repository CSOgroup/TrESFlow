/*
 * Module: DNA_QC
 * Runtime command:
 *   python3 bin/run_core_qc.py dna \
 *     --sb-group-map <sb_group_map.tsv> \
 *     --dna-mo-map <dna_mo_map.tsv> \
 *     --sample-counts <dna_sample_barcode.counts.tsv>... \
 *     --tag-records <dna_tag_records.tsv>... \
 *     --aligned-bam <sample_group_mark.bam>... \
 *     --nodup-bam <sample_group_mark_NoDup.bam>... \
 *     --outdir .
 *
 * Inputs:
 *   - derived sample-barcode group map
 *   - derived DNA modality map
 *   - published DNA sample-barcode counts tables
 *   - published DNA tag-record TSVs with SB / MO / CB annotations
 *   - published aligned DNA BAMs
 *   - published NoDup DNA BAMs
 * Outputs:
 *   - DNA QC stage-count tables
 *   - per-sample, per-group, and mark-aware DNA QC plots under dna/
 */

process DNA_QC {
    tag "dna_qc"

    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    path(sbGroupMap)
    path(dnaMoMap)
    path(sampleCounts)
    path(tagRecords)
    path(alignedBams)
    path(nodupBams)

    output:
    path("dna_sample_stage_counts.tsv"), emit: sample_table
    path("dna_group_stage_counts.tsv"), emit: group_table
    path("dna_group_mark_stage_counts.tsv"), emit: group_mark_table
    path("dna"), emit: plots

    script:
    def sampleCountArgs = sampleCounts.collect { "--sample-counts '${it}'" }.join(" \\\n      ")
    def tagRecordArgs = tagRecords.collect { "--tag-records '${it}'" }.join(" \\\n      ")
    def alignedBamArgs = alignedBams.collect { "--aligned-bam '${it}'" }.join(" \\\n      ")
    def nodupBamArgs = nodupBams.collect { "--nodup-bam '${it}'" }.join(" \\\n      ")

    """
    if [[ ! -x "\$SAMTOOLS_BIN" ]]; then
      echo "Missing configured DNA runtime executable: \$SAMTOOLS_BIN" >&2
      exit 1
    fi

    export MPLCONFIGDIR="\$(mktemp -d /tmp/mplconfig-dna-qc.XXXXXX)"
    trap 'rm -rf "\${MPLCONFIGDIR}"' EXIT

    "\$PYTHON3_BIN" "${projectDir}/bin/run_core_qc.py" dna \\
      --sb-group-map "${sbGroupMap}" \\
      --dna-mo-map "${dnaMoMap}" \\
      ${sampleCountArgs} \\
      ${tagRecordArgs} \\
      ${alignedBamArgs} \\
      ${nodupBamArgs} \\
      --outdir "."
    """
}
