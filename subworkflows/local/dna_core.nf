/*
 * Subworkflow: DNA_CORE
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw DNA I1 / I2 / R1 / R2 FASTQs
 *   - shared sample-barcode group map TSV used to derive the effective DNA SB whitelist
 *   - explicit DNA modality-barcode and ligation whitelists
 *   - DNA modality map TSV and the shared sample-barcode group map TSV for Split_ReadsV2 dna mode
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
    ch_sb_input = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, i2, r1, r2, sbGroupMap)
    }

    TAG_DNA_SAMPLE_BARCODE(ch_sb_input)

    ch_mo_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, i2, modalityWhitelist)
    }

    ch_mo_input = ch_mo_meta
        .join(TAG_DNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, i2, modalityWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, i2, taggedR1, taggedR2, modalityWhitelist)
        }

    TAG_DNA_MODALITY_BARCODE(ch_mo_input)

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
            def r1BySplit = WorkflowSupport.asPathList(splitR1s).collectEntries { path ->
                def splitName = path.getName().replaceFirst('_R1\\.fq\\.gz$', '')
                [(splitName): path]
            }
            def r2BySplit = WorkflowSupport.asPathList(splitR2s).collectEntries { path ->
                def splitName = path.getName().replaceFirst('_R2\\.fq\\.gz$', '')
                [(splitName): path]
            }

            def splitNames = (r1BySplit.keySet() + r2BySplit.keySet()).unique().sort()
            splitNames.collect { splitName ->
                if( !r1BySplit.containsKey(splitName) || !r2BySplit.containsKey(splitName) ) {
                    throw new IllegalStateException("Missing split FASTQ mate for DNA split '${splitName}'")
                }
                tuple(splitName, sampleId, meta, r1BySplit[splitName], r2BySplit[splitName])
            }
        }

    ch_align_rg = SPLIT_DNA_READS.out.rg_headers
        .flatMap { sampleId, meta, rgHeaders ->
            WorkflowSupport.asPathList(rgHeaders).collect { rgHeader ->
                def splitName = rgHeader.getName().replaceFirst('^SAM_RG_Header_', '').replaceFirst('\\.tsv$', '')
                tuple(splitName, sampleId, meta, rgHeader)
            }
        }

    ch_align_input = ch_align_fastqs
        .join(ch_align_rg)
        .map { splitName, sampleId, metaFromFastq, splitR1, splitR2, sampleIdFromRg, metaFromRg, rgHeader ->
            def suffix = splitName.replaceFirst("^${sampleId}_", '')
            def tokens = suffix.tokenize('_')
            if( tokens.size() < 2 ) {
                throw new IllegalStateException("Unable to derive DNA group and modality from split output '${splitName}'")
            }

            def group = tokens[0]
            def modality = tokens[1..-1].join('_')
            def sampleGroup = "${sampleId}_${group}"

            tuple(
                splitName,
                metaFromFastq,
                sampleGroup,
                modality,
                splitR1,
                splitR2,
                rgHeader,
                params.dna_bwa_reference as String,
                params.dna_blacklist_bed as String,
                (params.dna_effective_genome_size as String)
            )
        }

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
                (params.dna_effective_genome_size as String)
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
