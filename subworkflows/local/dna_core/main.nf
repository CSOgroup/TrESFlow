/*
 * Subworkflow: DNA_CORE
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw DNA I1 / I2 / R1 / R2 FASTQs
 *   - internally derived sample-barcode group map TSV used to derive the effective DNA SB whitelist
 *   - internally derived DNA modality-barcode whitelist plus the configured ligation whitelist
 *   - internally derived DNA modality map TSV and shared sample-barcode group map TSV for Split_ReadsV2 dna mode
 *   - DNA alignment resources carried through sample metadata
 * Outputs:
 *   - DNA FASTQs tagged with SB, MO, then CB comments
 *   - trim_galore paired-end FASTQs from the CB-tagged DNA reads
 *   - Split_ReadsV2 per-group per-mark DNA FASTQs and SAM RG headers
 *   - AlignDNA filtered BAMs, BAM indexes, and properly paired mapped-read counts
 *   - GATK duplicate-marked BAMs, BAM indexes, and duplicate metrics
 *   - duplicate-filtered NoDup BAMs and indexes
 *   - bigWig coverage tracks from the NoDup BAMs
 *   - barcode count/stat files from all wrapped DNA tagging steps
 */

import WorkflowSupport

include { TAG_DNA_SAMPLE_BARCODE }   from '../../modules/local/tag_dna_sb/main'
include { TAG_DNA_MODALITY_BARCODE } from '../../modules/local/tag_dna_modality/main'
include { TAG_DNA_CELL_BARCODE }     from '../../modules/local/tag_dna_cell_barcode/main'
include { TRIM_DNA_FASTQS }          from '../../modules/local/trim_dna_fastqs/main'
include { SPLIT_DNA_READS }          from '../../modules/local/split_dna_reads/main'
include { ALIGN_DNA }                from '../../modules/local/align_dna/main'
include { MARK_DUPLICATES_DNA }      from '../../modules/local/mark_duplicates_dna/main'
include { SPLIT_DUPLICATES_DNA }     from '../../modules/local/split_duplicates_dna/main'
include { BAM_COVERAGE_DNA }         from '../../modules/local/bam_coverage_dna/main'

workflow DNA_CORE {
    take:
    ch_dna_samples

    main:
    // Tag sample barcodes from the DNA I2 stream.
    ch_sb_input = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, i2, r1, r2, sbGroupMap)
    }

    TAG_DNA_SAMPLE_BARCODE(ch_sb_input)

    // Add modality barcodes from the same I2 read set.
    ch_mo_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, i2, modalityWhitelist)
    }

    ch_mo_input = ch_mo_meta
        .join(TAG_DNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, i2, modalityWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, i2, taggedR1, taggedR2, modalityWhitelist)
        }

    TAG_DNA_MODALITY_BARCODE(ch_mo_input)

    // Add ligation-derived cell barcodes from I1.
    ch_cb_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, i1, cellWhitelist)
    }

    ch_cb_input = ch_cb_meta
        .join(TAG_DNA_MODALITY_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, i1, cellWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, i1, taggedR1, taggedR2, cellWhitelist)
        }

    TAG_DNA_CELL_BARCODE(ch_cb_input)
    TRIM_DNA_FASTQS(TAG_DNA_CELL_BARCODE.out.tagged)

    // Split trimmed DNA reads by sample-barcode group and modality mark.
    ch_split_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, moMap, sbGroupMap)
    }

    ch_split_input = ch_split_meta
        .join(TRIM_DNA_FASTQS.out.trimmed)
        .map { sampleId, metaFromInput, moMap, sbGroupMap, metaFromTrim, trimmedR1, trimmedR2 ->
            tuple(sampleId, metaFromInput, trimmedR1, trimmedR2, moMap, sbGroupMap)
        }

    SPLIT_DNA_READS(ch_split_input)

    ch_align_fastqs = SPLIT_DNA_READS.out.split_fastqs
        .flatMap { sampleId, meta, splitR1s, splitR2s ->
            WorkflowSupport.pairDnaSplitFastqs(splitR1s, splitR2s).collect { split ->
                tuple(split.splitName, sampleId, meta, split.r1, split.r2)
            }
        }

    ch_align_rg = SPLIT_DNA_READS.out.rg_headers
        .flatMap { sampleId, meta, rgHeaders ->
            WorkflowSupport.collectDnaRgHeaders(rgHeaders).collect { rg ->
                tuple(rg.splitName, sampleId, meta, rg.rgHeader)
            }
        }

    ch_align_input = ch_align_fastqs
        .join(ch_align_rg)
        .map { splitName, sampleId, metaFromFastq, splitR1, splitR2, sampleIdFromRg, metaFromRg, rgHeader ->
            def splitMeta = WorkflowSupport.parseDnaSplitName(sampleId, splitName)

            tuple(
                splitName,
                metaFromFastq,
                splitMeta.sampleGroup,
                splitMeta.modality,
                splitR1,
                splitR2,
                rgHeader,
                metaFromFastq.dna_bwa_reference as String,
                metaFromFastq.dna_blacklist_bed as String,
                (metaFromFastq.dna_effective_genome_size as String)
            )
        }

    // Finish the DNA core with alignment, duplicate marking, NoDup extraction, and coverage.
    ALIGN_DNA(ch_align_input)
    MARK_DUPLICATES_DNA(ALIGN_DNA.out.bam)
    SPLIT_DUPLICATES_DNA(MARK_DUPLICATES_DNA.out.bam)

    ch_nodup_for_coverage = SPLIT_DUPLICATES_DNA.out.bam
        .join(SPLIT_DUPLICATES_DNA.out.bai)
        .map { splitName, metaFromBam, noDupBam, metaFromBai, noDupBai ->
            tuple(
                splitName,
                metaFromBam,
                noDupBam,
                noDupBai,
                (metaFromBam.dna_effective_genome_size as String)
            )
        }

    BAM_COVERAGE_DNA(ch_nodup_for_coverage)

    ch_barcode_reports = TAG_DNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_DNA_MODALITY_BARCODE.out.metrics)
        .mix(TAG_DNA_CELL_BARCODE.out.metrics)

    emit:
    tagged_fastqs   = TAG_DNA_CELL_BARCODE.out.tagged
    trimmed_fastqs  = TRIM_DNA_FASTQS.out.trimmed
    split_fastqs    = SPLIT_DNA_READS.out.split_fastqs
    rg_headers      = SPLIT_DNA_READS.out.rg_headers
    aligned_bams    = ALIGN_DNA.out.bam
    aligned_bais    = ALIGN_DNA.out.bai
    alignment_barcode_counts = ALIGN_DNA.out.barcode_counts
    markeddup_bams = MARK_DUPLICATES_DNA.out.bam
    markeddup_bais = MARK_DUPLICATES_DNA.out.bai
    duplicate_metrics = MARK_DUPLICATES_DNA.out.metrics
    nodup_bams = SPLIT_DUPLICATES_DNA.out.bam
    nodup_bais = SPLIT_DUPLICATES_DNA.out.bai
    coverage_bigwigs = BAM_COVERAGE_DNA.out.bw
    barcode_reports = ch_barcode_reports
}
