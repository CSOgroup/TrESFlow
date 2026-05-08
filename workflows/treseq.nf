/*
 * Workflow: TRESEQ
 * Core workflow:
 *   1. Parse a single hierarchical YAML samplesheet.
 *   2. Run the upstream RNA sample-barcode tagging step (Tag.codon) via a thin wrapper.
 *   3. Run the upstream RNA UMI tagging step (Tag_UMI.codon) via a thin wrapper.
 *   4. Run the upstream RNA cell-barcode tagging step (Tag_Lig3.codon) via a thin wrapper.
 *   5. Run the upstream RNA trim_galore step via a thin wrapper.
 *   6. Run the upstream RNA Split_ReadsV2 step in rna mode via a thin wrapper.
 *   7. Run the upstream RNA FqToSAM step via a thin wrapper.
 *   8. Run STARsolo alignment from grouped RNA unmapped SAMs.
 *   9. Generate filtered RNA BAMs from the STARsolo barcode calls.
 *  10. Generate stranded and unstranded RNA bigWigs from the filtered BAMs.
 *  11. Run the upstream DNA sample-barcode, modality-barcode, and cell-barcode tagging
 *      steps plus DNA trim_galore, Split_ReadsV2 dna mode, AlignDNA.sh,
 *      GATK MarkDuplicates, duplicate filtering to NoDup BAMs, and bamCoverage.
 */

import WorkflowSupport

include { RNA_CORE } from '../subworkflows/local/rna_core'
include { DNA_CORE } from '../subworkflows/local/dna_core'
include { SEQUENCING_EFFICIENCY } from '../modules/local/sequencing_efficiency/main'

def toRnaCoreInput(final Map row) {
    tuple(
        row.id,
        row,
        file(row.i1),
        file(row.r1),
        file(row.r2),
        file(row.cell_whitelist),
        file(row.sb_group_map)
    )
}

def toDnaCoreInput(final Map row) {
    tuple(
        row.id,
        row,
        file(row.i1),
        file(row.i2),
        file(row.r1),
        file(row.r2),
        file(row.modality_whitelist),
        file(row.cell_whitelist),
        file(row.mo_map),
        file(row.sb_group_map)
    )
}

def samplesheetParseOptions() {
    return [
        outdir                    : params.outdir ?: "${projectDir}/results",
        barcode_defaults          : params.barcode_defaults,
    ]
}

def uniqueFiles(final Collection paths) {
    return paths
        .findAll { it }
        .collect { file(it) }
        .unique { it.toString() }
}

def validateCoreResourceContract(final List<Map> rnaRows, final List<Map> dnaRows, final int maxCpus) {
    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }
}

workflow TRESEQ {
    take:
    sampleRows

    main:
    // Parse the single supported YAML contract into modality-specific work rows.
    final List<Map> rnaRows = sampleRows.findAll { row -> row.modality == 'rna' }
    final List<Map> dnaRows = sampleRows.findAll { row -> row.modality == 'dna' }
    final int maxCpus = params.max_cpus as int

    validateCoreResourceContract(rnaRows, dnaRows, maxCpus)

    Channel
        .fromList(rnaRows)
        .map { row -> toRnaCoreInput(row) }
        .set { ch_rna_samples }

    Channel
        .fromList(dnaRows)
        .map { row -> toDnaCoreInput(row) }
        .set { ch_dna_samples }

    // RNA and DNA run as independent branches under the same samplesheet contract.
    RNA_CORE(ch_rna_samples)
    DNA_CORE(ch_dna_samples)

    Channel
        .fromList(uniqueFiles(sampleRows.collect { row -> row.sb_group_map }))
        .collect()
        .ifEmpty([])
        .set { ch_efficiency_sb_group_maps }

    Channel
        .fromList(uniqueFiles(dnaRows.collect { row -> row.mo_map }))
        .collect()
        .ifEmpty([])
        .set { ch_efficiency_dna_mo_maps }

    SEQUENCING_EFFICIENCY(
        sampleRows ? sampleRows[0] : [:],
        RNA_CORE.out.tres_tag_records.collect().ifEmpty([]),
        RNA_CORE.out.aligned_filtered_bams.map { splitName, meta, bam -> bam }.collect().ifEmpty([]),
        DNA_CORE.out.tres_tag_records.collect().ifEmpty([]),
        DNA_CORE.out.markeddup_bams.map { splitName, meta, bam -> bam }.collect().ifEmpty([]),
        DNA_CORE.out.nodup_bams.map { splitName, meta, bam -> bam }.collect().ifEmpty([]),
        ch_efficiency_sb_group_maps,
        ch_efficiency_dna_mo_maps
    )

    emit:
    tagged_fastqs               = RNA_CORE.out.tagged_fastqs
    trimmed_fastqs              = RNA_CORE.out.trimmed_fastqs
    split_fastqs                = RNA_CORE.out.split_fastqs
    rg_headers                  = RNA_CORE.out.rg_headers
    usam_files                  = RNA_CORE.out.usam_files
    aligned_solo_dirs           = RNA_CORE.out.aligned_solo_dirs
    aligned_filtered_bams       = RNA_CORE.out.aligned_filtered_bams
    aligned_stranded_bigwigs    = RNA_CORE.out.aligned_stranded_bigwigs
    aligned_unstranded_bigwigs = RNA_CORE.out.aligned_unstranded_bigwigs
    barcode_reports             = RNA_CORE.out.barcode_reports
    dna_tagged_fastqs           = DNA_CORE.out.tagged_fastqs
    dna_trimmed_fastqs          = DNA_CORE.out.trimmed_fastqs
    dna_split_fastqs            = DNA_CORE.out.split_fastqs
    dna_rg_headers              = DNA_CORE.out.rg_headers
    dna_aligned_bams            = DNA_CORE.out.aligned_bams
    dna_aligned_bais            = DNA_CORE.out.aligned_bais
    dna_alignment_barcode_counts = DNA_CORE.out.alignment_barcode_counts
    dna_markeddup_bams          = DNA_CORE.out.markeddup_bams
    dna_markeddup_bais          = DNA_CORE.out.markeddup_bais
    dna_duplicate_metrics       = DNA_CORE.out.duplicate_metrics
    dna_nodup_bams              = DNA_CORE.out.nodup_bams
    dna_nodup_bais              = DNA_CORE.out.nodup_bais
    dna_coverage_bigwigs        = DNA_CORE.out.coverage_bigwigs
    dna_coverage_warnings       = DNA_CORE.out.coverage_warnings
    dna_barcode_reports         = DNA_CORE.out.barcode_reports
    sequencing_efficiency_reports = SEQUENCING_EFFICIENCY.out.reports
    sequencing_efficiency_warnings = SEQUENCING_EFFICIENCY.out.warnings
}
