/*
 * Subworkflow: INITIAL_RNA_TAGGING
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw RNA R1 / R2 FASTQs
 *   - RNA sample-barcode whitelist
 * Outputs:
 *   - RNA FASTQs tagged with SB then UM comments
 *   - barcode count/stat files from both wrapped upstream steps
 */

include { TAG_RNA_SAMPLE_BARCODE } from '../../modules/local/tag_rna_sb/main'
include { TAG_RNA_UMI }            from '../../modules/local/tag_rna_umi/main'

workflow INITIAL_RNA_TAGGING {
    take:
    ch_rna_samples

    main:
    TAG_RNA_SAMPLE_BARCODE(ch_rna_samples)

    ch_raw_r2 = ch_rna_samples.map { sampleId, meta, r1, r2, whitelist ->
        tuple(sampleId, meta, r2)
    }

    ch_umi_input = ch_raw_r2
        .join(TAG_RNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, rawR2, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, rawR2, taggedR1, taggedR2)
        }

    TAG_RNA_UMI(ch_umi_input)

    ch_barcode_reports = TAG_RNA_SAMPLE_BARCODE.out.metrics.mix(TAG_RNA_UMI.out.metrics)

    emit:
    tagged_fastqs   = TAG_RNA_UMI.out.tagged
    barcode_reports = ch_barcode_reports
}
