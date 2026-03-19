/*
 * Subworkflow: INITIAL_RNA_TAGGING
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw RNA I1 / R1 / R2 FASTQs
 *   - RNA sample-barcode and cell-barcode whitelists
 * Outputs:
 *   - RNA FASTQs tagged with SB, UM, then CB comments
 *   - trim_galore paired-end FASTQs from the CB-tagged reads
 *   - barcode count/stat files from all wrapped upstream steps
 */

include { TAG_RNA_SAMPLE_BARCODE } from '../../modules/local/tag_rna_sb/main'
include { TAG_RNA_UMI }            from '../../modules/local/tag_rna_umi/main'
include { TAG_RNA_CELL_BARCODE }   from '../../modules/local/tag_rna_cell_barcode/main'
include { TRIM_RNA_FASTQS }        from '../../modules/local/trim_rna_fastqs/main'

workflow INITIAL_RNA_TAGGING {
    take:
    ch_rna_samples

    main:
    ch_sb_input = ch_rna_samples.map { sampleId, meta, i1, r1, r2, sampleWhitelist, cellWhitelist ->
        tuple(sampleId, meta, r1, r2, sampleWhitelist)
    }

    TAG_RNA_SAMPLE_BARCODE(ch_sb_input)

    ch_raw_r2 = ch_rna_samples.map { sampleId, meta, i1, r1, r2, sampleWhitelist, cellWhitelist ->
        tuple(sampleId, meta, r2)
    }

    ch_umi_input = ch_raw_r2
        .join(TAG_RNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, rawR2, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, rawR2, taggedR1, taggedR2)
        }

    TAG_RNA_UMI(ch_umi_input)

    ch_cb_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, sampleWhitelist, cellWhitelist ->
        tuple(sampleId, meta, i1, cellWhitelist)
    }

    ch_cb_input = ch_cb_meta
        .join(TAG_RNA_UMI.out.tagged)
        .map { sampleId, metaFromInput, i1, cellWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, i1, taggedR1, taggedR2, cellWhitelist)
        }

    TAG_RNA_CELL_BARCODE(ch_cb_input)
    TRIM_RNA_FASTQS(TAG_RNA_CELL_BARCODE.out.tagged)

    ch_barcode_reports = TAG_RNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_RNA_UMI.out.metrics)
        .mix(TAG_RNA_CELL_BARCODE.out.metrics)

    emit:
    tagged_fastqs    = TAG_RNA_CELL_BARCODE.out.tagged
    trimmed_fastqs   = TRIM_RNA_FASTQS.out.trimmed
    barcode_reports  = ch_barcode_reports
}
