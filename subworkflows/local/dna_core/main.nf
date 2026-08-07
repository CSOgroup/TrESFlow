/*
 * Subworkflow: DNA_CORE
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw DNA I1 / I2 / R1 / R2 FASTQs
 *   - internally derived DNA sample-barcode group map TSV used to derive the effective DNA SB whitelist
 *   - internally derived DNA modality-barcode whitelist plus the configured ligation whitelist
 *   - internally derived DNA modality map TSV and shared sample-barcode group map TSV for Split_ReadsV2 dna mode
 *   - DNA alignment resources carried through sample metadata
 * Outputs:
 *   - DNA FASTQs tagged with SB, MO, then CB comments
 *   - uncompressed trim_galore paired-end FASTQs from the CB-tagged DNA reads
 *   - Split_ReadsV2 per-group per-mark DNA FASTQs and SAM RG headers
 *   - AlignDNA filtered BAMs and BAM indexes
 *   - GATK duplicate-marked BAMs, BAM indexes, and duplicate metrics
 *   - duplicate-filtered NoDup BAMs and indexes
 *   - bigWig coverage tracks from the NoDup BAMs
 *   - barcode count/stat files from all wrapped DNA tagging steps
 */

import WorkflowSupport

include { TAG_DNA_SAMPLE_BARCODE }   from '../../../modules/local/tag_dna_sb/main'
include { TAG_DNA_MODALITY_BARCODE } from '../../../modules/local/tag_dna_modality/main'
include { TAG_DNA_CELL_BARCODE }     from '../../../modules/local/tag_dna_cell_barcode/main'
include { TRIM_DNA_FASTQS }          from '../../../modules/local/trim_dna_fastqs/main'
include { SPLIT_DNA_READS }          from '../../../modules/local/split_dna_reads/main'
include { ALIGN_DNA }                from '../../../modules/local/align_dna/main'
include { GATK4_MARKDUPLICATES }     from '../../../modules/nf-core/gatk4/markduplicates/main'
include { NORMALIZE_DNA_MARKDUPLICATES } from '../../../modules/local/normalize_dna_markduplicates/main'
include { SPLIT_DUPLICATES_DNA }     from '../../../modules/local/split_duplicates_dna/main'
include { CHECK_DNA_NODUP_BAM }      from '../../../modules/local/check_dna_nodup_bam/main'
include { DEEPTOOLS_BAMCOVERAGE }    from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { NORMALIZE_DNA_BAMCOVERAGE } from '../../../modules/local/normalize_dna_bamcoverage/main'

def selectDnaIndexRead(final Map meta, final Object i1, final Object i2, final String fieldName) {
    return meta[fieldName] == 'i1' ? i1 : i2
}

def selectDnaLigationRead(final Map meta, final Object i1, final Object i2) {
    return i1
}

def dnaLigationStartPositions(final Map meta) {
    return meta.dna_tagmentation == 'dual' ? '41,79,117' : '15,53,91'
}

def nfcoreDnaMeta(final String splitName, final Map meta, final String stage) {
    return meta + [
        id             : splitName,
        tres_sample_id : meta.id,
        tres_split_name: splitName,
        tres_stage     : stage,
    ]
}

def restoreDnaMeta(final Map meta) {
    return meta + [id: meta.tres_sample_id ?: meta.id]
}

workflow DNA_CORE {
    take:
    ch_dna_samples

    main:
    ch_versions = channel.empty()

    // Tag sample barcodes from the tagmentation-specific DNA index stream.
    ch_sb_input = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        def sbIndexRead = selectDnaIndexRead(meta, i1, i2, 'dna_sample_index_read')
        tuple(sampleId, meta, sbIndexRead, r1, r2, sbGroupMap)
    }

    TAG_DNA_SAMPLE_BARCODE(ch_sb_input)
    ch_versions = ch_versions.mix(TAG_DNA_SAMPLE_BARCODE.out.versions)

    // Add modality barcodes from the tagmentation-specific DNA index stream.
    ch_mo_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        def moIndexRead = selectDnaIndexRead(meta, i1, i2, 'dna_modality_index_read')
        tuple(sampleId, meta, moIndexRead, modalityWhitelist)
    }

    ch_mo_input = ch_mo_meta
        .join(TAG_DNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, indexRead, modalityWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, indexRead, taggedR1, taggedR2, modalityWhitelist)
        }

    TAG_DNA_MODALITY_BARCODE(ch_mo_input)
    ch_versions = ch_versions.mix(TAG_DNA_MODALITY_BARCODE.out.versions)

    // Add ligation-derived cell barcodes from DNA I1. Legacy single-tagmentation
    // data still keeps SB/MO on I2, but L1/L2/L3 are read from I1.
    ch_cb_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        def ligationRead = selectDnaLigationRead(meta, i1, i2)
        def ligationMeta = meta + [
            dna_ligation_index_read    : 'i1',
            dna_ligation_start_positions: dnaLigationStartPositions(meta),
        ]
        tuple(sampleId, ligationMeta, ligationRead, cellWhitelist)
    }

    ch_cb_input = ch_cb_meta
        .join(TAG_DNA_MODALITY_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, ligationRead, cellWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, ligationRead, taggedR1, taggedR2, cellWhitelist)
        }

    TAG_DNA_CELL_BARCODE(ch_cb_input)
    ch_versions = ch_versions.mix(TAG_DNA_CELL_BARCODE.out.versions)
    TRIM_DNA_FASTQS(TAG_DNA_CELL_BARCODE.out.tagged)
    ch_versions = ch_versions.mix(TRIM_DNA_FASTQS.out.versions)

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
    ch_versions = ch_versions.mix(SPLIT_DNA_READS.out.versions)

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
    ch_versions = ch_versions.mix(ALIGN_DNA.out.versions)

    ch_gatk_markduplicates_input = ALIGN_DNA.out.bam.map { splitName, meta, alignedBam ->
        tuple(nfcoreDnaMeta(splitName, meta, 'markeddup'), alignedBam)
    }

    GATK4_MARKDUPLICATES(ch_gatk_markduplicates_input, [], [])
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions_gatk4)
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions_samtools)

    ch_normalize_markduplicates_input = GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai)
        .join(GATK4_MARKDUPLICATES.out.metrics)
        .map { nfMeta, markedDupBam, markedDupBai, markedDupMetrics ->
            tuple(
                nfMeta.tres_split_name as String,
                restoreDnaMeta(nfMeta),
                markedDupBam,
                markedDupBai,
                markedDupMetrics
            )
        }

    NORMALIZE_DNA_MARKDUPLICATES(ch_normalize_markduplicates_input)
    ch_versions = ch_versions.mix(NORMALIZE_DNA_MARKDUPLICATES.out.versions)

    SPLIT_DUPLICATES_DNA(NORMALIZE_DNA_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(SPLIT_DUPLICATES_DNA.out.versions)

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

    CHECK_DNA_NODUP_BAM(ch_nodup_for_coverage)
    ch_versions = ch_versions.mix(CHECK_DNA_NODUP_BAM.out.versions)

    ch_deeptools_bamcoverage_input = CHECK_DNA_NODUP_BAM.out.ready.map { splitName, meta, noDupBam, noDupBai, effectiveGenomeSize ->
        tuple(nfcoreDnaMeta(splitName, meta, 'nodup_coverage'), noDupBam, noDupBai)
    }

    DEEPTOOLS_BAMCOVERAGE(
        ch_deeptools_bamcoverage_input,
        [],
        [],
        tuple([id: 'no_blacklist'], [])
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions_deeptools)
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions_samtools)

    ch_normalize_bamcoverage_input = DEEPTOOLS_BAMCOVERAGE.out.bigwig.map { nfMeta, bigwig ->
        tuple(
            nfMeta.tres_split_name as String,
            restoreDnaMeta(nfMeta),
            bigwig
        )
    }

    NORMALIZE_DNA_BAMCOVERAGE(ch_normalize_bamcoverage_input)
    ch_versions = ch_versions.mix(NORMALIZE_DNA_BAMCOVERAGE.out.versions)

    ch_barcode_reports = TAG_DNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_DNA_MODALITY_BARCODE.out.metrics)
        .mix(TAG_DNA_CELL_BARCODE.out.metrics)

    ch_barcode_report_files = TAG_DNA_SAMPLE_BARCODE.out.metrics
        .flatMap { sampleId, counts, stats -> [counts, stats] }
        .mix(TAG_DNA_MODALITY_BARCODE.out.metrics.flatMap { sampleId, counts, stats -> [counts, stats] })
        .mix(TAG_DNA_CELL_BARCODE.out.metrics.flatMap { sampleId, counts, tagRecords, statsL1, statsL2, statsL3 ->
            [counts, statsL1, statsL2, statsL3]
        })

    emit:
    tagged_fastqs   = TAG_DNA_CELL_BARCODE.out.tagged
    trimmed_fastqs  = TRIM_DNA_FASTQS.out.trimmed
    split_fastqs    = SPLIT_DNA_READS.out.split_fastqs
    rg_headers      = SPLIT_DNA_READS.out.rg_headers
    aligned_bams    = ALIGN_DNA.out.bam
    aligned_bais    = ALIGN_DNA.out.bai
    markeddup_bams = NORMALIZE_DNA_MARKDUPLICATES.out.bam
    markeddup_bais = NORMALIZE_DNA_MARKDUPLICATES.out.bai
    duplicate_metrics = NORMALIZE_DNA_MARKDUPLICATES.out.metrics
    nodup_bams = SPLIT_DUPLICATES_DNA.out.bam
    nodup_bais = SPLIT_DUPLICATES_DNA.out.bai
    coverage_bigwigs = NORMALIZE_DNA_BAMCOVERAGE.out.bw
    coverage_warnings = CHECK_DNA_NODUP_BAM.out.warnings
    barcode_reports = ch_barcode_reports
    barcode_report_files = ch_barcode_report_files
    tres_tag_records = TAG_DNA_CELL_BARCODE.out.tres_tag_records
    versions = ch_versions
}
