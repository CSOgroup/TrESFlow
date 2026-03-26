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
 *  12. Summarize the completed RNA and DNA core outputs into user-facing QC tables and plots.
 */

import WorkflowSupport

include { RNA_CORE } from '../subworkflows/local/rna_core'
include { DNA_CORE } from '../subworkflows/local/dna_core'
include { CORE_QC }  from '../subworkflows/local/core_qc'

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
        outdir                    : params.outdir,
        ligation_barcode_whitelist: params.ligation_barcode_whitelist,
        barcode_defaults          : params.barcode_defaults,
        rna_ref_base_dir          : params.rna_ref_base_dir,
        rna_align_species         : params.rna_align_species,
        dna_bwa_reference         : params.dna_bwa_reference,
        dna_blacklist_bed         : params.dna_blacklist_bed,
        dna_effective_genome_size : params.dna_effective_genome_size,
    ]
}

def validateCoreResourceContract(final List<Map> rnaRows, final List<Map> dnaRows, final int maxCpus) {
    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }

    if( rnaRows ) {
        WorkflowSupport.validateRnaAlignment(
            rnaRows[0].rna_ref_base_dir as String,
            rnaRows[0].rna_align_species as String
        )
    }

    if( dnaRows ) {
        try {
            WorkflowSupport.validateDnaAlignment(
                dnaRows[0].dna_bwa_reference as String,
                dnaRows[0].dna_blacklist_bed as String,
                dnaRows[0].dna_effective_genome_size as String
            )
        }
        catch( IllegalArgumentException e ) {
            error e.message
        }
    }
}

workflow TRESEQ {
    main:
    if( !params.samplesheet ) {
        error "Missing required parameter: --samplesheet"
    }

    // Parse the single supported YAML contract into modality-specific work rows.
    final List<Map> sampleRows = SamplesheetParser.parse(params.samplesheet as String, samplesheetParseOptions())
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

    final def ch_sb_group_map = Channel.value(file(sampleRows[0].sb_group_map))
    final def ch_dna_mo_map = dnaRows ? Channel.value(file(dnaRows[0].mo_map)) : Channel.empty()

    final def ch_rna_sample_counts = rnaRows
        ? RNA_CORE.out.sample_barcode_counts.map { sampleId, counts -> counts }.collect()
        : Channel.empty()
    final def ch_rna_cell_counts = rnaRows
        ? RNA_CORE.out.cell_barcode_counts.map { sampleId, counts -> counts }.collect()
        : Channel.empty()
    final def ch_rna_solo_dirs = rnaRows
        ? RNA_CORE.out.aligned_solo_dirs.map { splitName, meta, soloDir -> soloDir }.collect()
        : Channel.empty()

    final def ch_dna_sample_counts = dnaRows
        ? DNA_CORE.out.sample_barcode_counts.map { sampleId, counts -> counts }.collect()
        : Channel.empty()
    final def ch_dna_tag_records = dnaRows
        ? DNA_CORE.out.cell_barcode_tag_records.map { sampleId, tagRecords -> tagRecords }.collect()
        : Channel.empty()
    final def ch_dna_aligned_bams = dnaRows
        ? DNA_CORE.out.aligned_bams.map { splitName, meta, bam -> bam }.collect()
        : Channel.empty()
    final def ch_dna_nodup_bams = dnaRows
        ? DNA_CORE.out.nodup_bams.map { splitName, meta, bam -> bam }.collect()
        : Channel.empty()

    CORE_QC(
        ch_sb_group_map,
        ch_rna_sample_counts,
        ch_rna_cell_counts,
        ch_rna_solo_dirs,
        ch_dna_mo_map,
        ch_dna_sample_counts,
        ch_dna_tag_records,
        ch_dna_aligned_bams,
        ch_dna_nodup_bams
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
    aligned_unstranded_bigwigs  = RNA_CORE.out.aligned_unstranded_bigwigs
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
    dna_barcode_reports         = DNA_CORE.out.barcode_reports
    qc_rna_sample_table         = CORE_QC.out.rna_sample_table
    qc_rna_group_table          = CORE_QC.out.rna_group_table
    qc_rna_plots                = CORE_QC.out.rna_plots
    qc_dna_sample_table         = CORE_QC.out.dna_sample_table
    qc_dna_group_table          = CORE_QC.out.dna_group_table
    qc_dna_group_mark_table     = CORE_QC.out.dna_group_mark_table
    qc_dna_plots                = CORE_QC.out.dna_plots
}
