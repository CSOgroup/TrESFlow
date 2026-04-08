/*
 * Subworkflow: RNA_CORE
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw RNA I1 / R1 / R2 FASTQs
 *   - RNA cell-barcode whitelist
 *   - internally derived sample-barcode group map TSV, used both to derive the
 *     effective RNA sample-barcode whitelist and to split grouped RNA reads
 *   - RNA alignment reference base dir carried through sample metadata
 * Outputs:
 *   - RNA FASTQs tagged with SB, UM, then CB comments
 *   - trim_galore paired-end FASTQs from the CB-tagged reads
 *   - Split_ReadsV2 per-group RNA FASTQs and SAM RG headers
 *   - FqToSAM unmapped SAM files from each split RNA FASTQ pair
 *   - STARsolo outputs, filtered BAMs, and bigWigs from the decomposed RNA alignment path
 *   - barcode count/stat files from all wrapped RNA steps
 */

import WorkflowSupport

include { TAG_RNA_SAMPLE_BARCODE } from '../../../modules/local/tag_rna_sb/main'
include { TAG_RNA_UMI }            from '../../../modules/local/tag_rna_umi/main'
include { TAG_RNA_CELL_BARCODE }   from '../../../modules/local/tag_rna_cell_barcode/main'
include { TRIM_RNA_FASTQS }        from '../../../modules/local/trim_rna_fastqs/main'
include { SPLIT_RNA_READS }        from '../../../modules/local/split_rna_reads/main'
include { FQ_TO_SAM }              from '../../../modules/local/fq_to_sam/main'
include { RNA_STARSOLO_ALIGN }     from '../../../modules/local/rna_starsolo_align/main'
include { RNA_FILTERED_BAM }       from '../../../modules/local/rna_filtered_bam/main'
include { RNA_COVERAGE }           from '../../../modules/local/rna_coverage/main'

workflow RNA_CORE {
    take:
    ch_rna_samples

    main:
    ch_versions = channel.empty()

    // Tag sample barcodes from the RNA read-2 stream.
    ch_sb_input = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, r1, r2, sbGroupMap)
    }

    TAG_RNA_SAMPLE_BARCODE(ch_sb_input)
    ch_versions = ch_versions.mix(TAG_RNA_SAMPLE_BARCODE.out.versions)

    // Add UMIs after sample-barcode tagging.
    ch_raw_r2 = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, r2)
    }

    ch_umi_input = ch_raw_r2
        .join(TAG_RNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, rawR2, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, rawR2, taggedR1, taggedR2)
        }

    TAG_RNA_UMI(ch_umi_input)
    ch_versions = ch_versions.mix(TAG_RNA_UMI.out.versions)

    // Add ligation-derived cell barcodes from I1.
    ch_cb_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, i1, cellWhitelist)
    }

    ch_cb_input = ch_cb_meta
        .join(TAG_RNA_UMI.out.tagged)
        .map { sampleId, metaFromInput, i1, cellWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, i1, taggedR1, taggedR2, cellWhitelist)
        }

    TAG_RNA_CELL_BARCODE(ch_cb_input)
    ch_versions = ch_versions.mix(TAG_RNA_CELL_BARCODE.out.versions)
    TRIM_RNA_FASTQS(TAG_RNA_CELL_BARCODE.out.tagged)
    ch_versions = ch_versions.mix(TRIM_RNA_FASTQS.out.versions)

    // Split trimmed reads by sample-barcode group before FQ_TO_SAM.
    ch_split_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, sbGroupMap)
    }

    ch_split_input = ch_split_meta
        .join(TRIM_RNA_FASTQS.out.trimmed)
        .map { sampleId, metaFromInput, sbGroupMap, metaFromTrim, trimmedR1, trimmedR2 ->
            tuple(sampleId, metaFromInput, trimmedR1, trimmedR2, sbGroupMap)
        }

    SPLIT_RNA_READS(ch_split_input)
    ch_versions = ch_versions.mix(SPLIT_RNA_READS.out.versions)

    ch_fq_to_sam_input = SPLIT_RNA_READS.out.split_fastqs
        .flatMap { sampleId, meta, splitR1s, splitR2s ->
            WorkflowSupport.pairRnaSplitFastqs(sampleId, splitR1s, splitR2s).collect { split ->
                tuple(split.splitName, meta, split.r1, split.r2)
            }
        }

    FQ_TO_SAM(ch_fq_to_sam_input)
    ch_versions = ch_versions.mix(FQ_TO_SAM.out.versions)

    // Align each grouped unmapped SAM independently with STARsolo.
    ch_starsolo_input = FQ_TO_SAM.out.usam
        .map { splitName, meta, usam ->
            tuple(
                splitName,
                meta,
                usam,
                meta.rna_ref_base_dir as String,
                meta.rna_align_species as String
            )
        }

    RNA_STARSOLO_ALIGN(ch_starsolo_input)
    ch_versions = ch_versions.mix(RNA_STARSOLO_ALIGN.out.versions)

    ch_filtered_bam_input = RNA_STARSOLO_ALIGN.out.solo_dir
        .join(RNA_STARSOLO_ALIGN.out.aligned_bam)
        .map { splitName, metaFromSolo, soloDir, metaFromBam, alignedBam ->
            tuple(splitName, metaFromSolo, soloDir, alignedBam)
        }

    RNA_FILTERED_BAM(ch_filtered_bam_input)
    ch_versions = ch_versions.mix(RNA_FILTERED_BAM.out.versions)

    ch_coverage_input = RNA_FILTERED_BAM.out.filtered_bam
        .map { splitName, meta, filteredBam ->
            tuple(
                splitName,
                meta,
                filteredBam,
                meta.rna_ref_base_dir as String,
                meta.rna_align_species as String
            )
        }

    RNA_COVERAGE(ch_coverage_input)
    ch_versions = ch_versions.mix(RNA_COVERAGE.out.versions)

    ch_barcode_reports = TAG_RNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_RNA_UMI.out.metrics)
        .mix(TAG_RNA_CELL_BARCODE.out.metrics)

    emit:
    tagged_fastqs    = TAG_RNA_CELL_BARCODE.out.tagged
    trimmed_fastqs   = TRIM_RNA_FASTQS.out.trimmed
    split_fastqs     = SPLIT_RNA_READS.out.split_fastqs
    rg_headers       = SPLIT_RNA_READS.out.rg_headers
    usam_files       = FQ_TO_SAM.out.usam
    aligned_solo_dirs = RNA_STARSOLO_ALIGN.out.solo_dir
    aligned_filtered_bams = RNA_FILTERED_BAM.out.filtered_bam
    aligned_stranded_bigwigs = RNA_COVERAGE.out.stranded_bw
    aligned_unstranded_bigwigs = RNA_COVERAGE.out.unstranded_bw
    barcode_reports  = ch_barcode_reports
    versions         = ch_versions
}
