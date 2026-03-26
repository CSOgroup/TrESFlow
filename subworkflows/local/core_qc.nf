/*
 * Subworkflow: CORE_QC
 * Inputs:
 *   - derived sample-barcode group map
 *   - derived DNA modality map when DNA is present
 *   - collected RNA QC inputs from the current RNA core outputs
 *   - collected DNA QC inputs from the current DNA core outputs
 * Outputs:
 *   - RNA and DNA QC tables plus plot directories under qc/
 */

include { RNA_QC } from '../../modules/local/rna_qc/main'
include { DNA_QC } from '../../modules/local/dna_qc/main'

workflow CORE_QC {
    take:
    ch_sb_group_map
    ch_rna_sample_counts
    ch_rna_cell_counts
    ch_rna_solo_dirs
    ch_dna_mo_map
    ch_dna_sample_counts
    ch_dna_tag_records
    ch_dna_aligned_bams
    ch_dna_nodup_bams

    main:
    RNA_QC(ch_sb_group_map, ch_rna_sample_counts, ch_rna_cell_counts, ch_rna_solo_dirs)
    DNA_QC(ch_sb_group_map, ch_dna_mo_map, ch_dna_sample_counts, ch_dna_tag_records, ch_dna_aligned_bams, ch_dna_nodup_bams)

    emit:
    rna_sample_table = RNA_QC.out.sample_table
    rna_group_table = RNA_QC.out.group_table
    rna_plots = RNA_QC.out.plots
    dna_sample_table = DNA_QC.out.sample_table
    dna_group_table = DNA_QC.out.group_table
    dna_group_mark_table = DNA_QC.out.group_mark_table
    dna_plots = DNA_QC.out.plots
}
