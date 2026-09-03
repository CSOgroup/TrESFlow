/*
 * Subworkflow: DNA_CORE
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - ordered raw DNA I1 / I2 / R1 / R2 FASTQ collections (one logical library row)
 *   - internally derived DNA sample-barcode group map TSV used to derive the effective DNA SB whitelist
 *   - internally derived DNA modality-barcode whitelist plus the configured ligation whitelist
 *   - internally derived DNA modality map TSV and shared sample-barcode group map TSV for Split_ReadsV2 dna mode
 *   - DNA alignment resources carried through sample metadata
 * Outputs:
 *   - DNA FASTQs tagged with SB, MO, then CB comments
 *   - uncompressed trim_galore paired-end FASTQs from the CB-tagged DNA reads
 *   - plain Split_ReadsV2 per-group per-mark DNA FASTQs for computation plus optional independently compressed publication copies and SAM RG headers
 *   - AlignDNA filtered BAMs and BAM indexes
 *   - GATK duplicate-marked BAMs, BAM indexes, and duplicate metrics
 *   - duplicate-filtered NoDup BAMs and indexes
 *   - bigWig coverage tracks from the NoDup BAMs
 *   - barcode count/stat files from all wrapped DNA tagging steps
 */

include { TAG_DNA_SAMPLE_BARCODE }   from '../../../modules/local/tag_dna_sb/main'
include { TAG_DNA_MODALITY_BARCODE } from '../../../modules/local/tag_dna_modality/main'
include { TAG_DNA_CELL_BARCODE }     from '../../../modules/local/tag_dna_cell_barcode/main'
include { DUAL_TAG_ARTIFACT_FILTER } from '../../../modules/local/dual_tag_artifact_filter/main'
include { TRIM_DNA_FASTQS }          from '../../../modules/local/trim_dna_fastqs/main'
include { BARCODE_GATE_METRICS }      from '../../../modules/local/barcode_gate_metrics/main'
include { SPLIT_DNA_READS }          from '../../../modules/local/split_dna_reads/main'
include { COMPRESS_SPLIT_FASTQS as COMPRESS_DNA_SPLIT_FASTQS } from '../../../modules/local/compress_split_fastqs/main'
include { ALIGN_DNA }                from '../../../modules/local/align_dna/main'
include { FILTER_CANONICAL_DNA_ALIGNED_BAM } from '../../../modules/local/filter_canonical_dna_aligned_bam/main'
include { GATK4_MARKDUPLICATES }     from '../../../modules/nf-core/gatk4/markduplicates/main'
include { NORMALIZE_DNA_MARKDUPLICATES } from '../../../modules/local/normalize_dna_markduplicates/main'
include { SPLIT_DUPLICATES_DNA }     from '../../../modules/local/split_duplicates_dna/main'
include { DEEPTOOLS_BAMCOVERAGE }    from '../../../modules/nf-core/deeptools/bamcoverage/main'

def asDnaPathList(value) {
    return value instanceof List ? value : [value]
}

def pairDnaSplitFastqs(splitR1s, splitR2s) {
    def r1BySplit = asDnaPathList(splitR1s).collectEntries { path ->
        [(path.getName().replaceFirst('_R1\\.(?:fastq|fq)(?:\\.gz)?$', '')): path]
    }
    def r2BySplit = asDnaPathList(splitR2s).collectEntries { path ->
        [(path.getName().replaceFirst('_R2\\.(?:fastq|fq)(?:\\.gz)?$', '')): path]
    }

    def splitNames = (r1BySplit.keySet() + r2BySplit.keySet()).unique().sort()
    splitNames.collect { splitName ->
        if( !r1BySplit.containsKey(splitName) || !r2BySplit.containsKey(splitName) ) {
            throw new IllegalStateException("Missing split FASTQ mate for DNA split '${splitName}'")
        }
        [splitName: splitName, r1: r1BySplit[splitName], r2: r2BySplit[splitName]]
    }
}

def collectDnaRgHeaders(rgHeaders) {
    asDnaPathList(rgHeaders).collect { rgHeader ->
        [
            splitName: rgHeader.getName()
                .replaceFirst('^SAM_RG_Header_', '')
                .replaceFirst('\\.tsv$', ''),
            rgHeader : rgHeader,
        ]
    }
}

def parseDnaSplitName(sampleId, splitName) {
    def suffix = splitName.replaceFirst("^${sampleId}_", '')
    def tokens = suffix.tokenize('_')
    if( tokens.size() < 2 ) {
        throw new IllegalStateException(
            "Unable to derive DNA group and modality from split output '${splitName}'"
        )
    }

    def group = tokens[0]
    [
        group      : group,
        modality   : tokens[1..-1].join('_'),
        sampleGroup: "${sampleId}_${group}",
    ]
}

def selectDnaIndexRead(meta, i1, i2, fieldName) {
    return meta[fieldName] == 'i1' ? i1 : i2
}

def selectDnaLigationRead(meta, i1, i2) {
    return i1
}

def dnaLigationStartPositions(meta) {
    return meta.dna_tagmentation == 'dual' ? '41,79,117' : '15,53,91'
}

def nfcoreDnaMeta(splitName, meta, stage) {
    return meta + [
        id             : splitName,
        tres_sample_id : meta.id,
        tres_split_name: splitName,
        tres_stage     : stage,
    ]
}

def restoreDnaMeta(meta) {
    return meta + [id: meta.tres_sample_id ?: meta.id]
}

def hasMappedDnaReads(mappedReadsFile) {
    def value = mappedReadsFile.text.trim()
    if( !(value ==~ /[0-9]+/) ) {
        throw new IllegalStateException(
            "Invalid NoDup mapped-read count '${value}' in ${mappedReadsFile}"
        )
    }
    return value.toLong() > 0
}

workflow DNA_CORE {
    take:
    ch_dna_samples

    main:
    ch_versions = channel.empty()
    def codonScriptsRoot = java.lang.System.getProperty('tresflow.resolvedCoreScriptsDir')

    def tagHelperScripts = [
        file("${projectDir}/bin/run_tag.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def tagCodonScripts = [
        file("${codonScriptsRoot}/Tag.codon", checkIfExists: true),
        file("${codonScriptsRoot}/utils.codon", checkIfExists: true),
    ]
    def cellHelperScripts = [
        file("${projectDir}/bin/run_tag_lig3.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def cellCodonScripts = [
        file("${codonScriptsRoot}/Tag_Lig3.codon", checkIfExists: true),
        file("${codonScriptsRoot}/utils.codon", checkIfExists: true),
    ]
    def splitHelperScripts = [
        file("${projectDir}/bin/run_split_reads_dna.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def splitCodonScripts = [
        file("${codonScriptsRoot}/Split_ReadsV2.codon", checkIfExists: true),
    ]
    def trimHelperScript = file("${projectDir}/bin/run_trim_galore.py", checkIfExists: true)
    def dualTagFilterHelperScript = file(
        "${projectDir}/bin/run_dual_tag_artifact_filter.py",
        checkIfExists: true
    )
    def alignRuntimeScripts = [
        file("${projectDir}/scripts/core_runtime/AlignDNA.sh", checkIfExists: true),
    ]
    def canonicalBamRuntimeScripts = [
        file("${projectDir}/scripts/core_runtime/FilterCanonicalBam.sh", checkIfExists: true),
    ]
    def splitDuplicatesRuntimeScripts = [
        file("${projectDir}/scripts/core_runtime/SplitDuplicatesDNA.sh", checkIfExists: true),
    ]
    def bwaMem2IndexSuffixes = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']

    // Tag sample barcodes from the tagmentation-specific DNA index stream.
    ch_sb_input = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        def sbIndexRead = selectDnaIndexRead(meta, i1, i2, 'dna_sample_index_read')
        tuple(sampleId, meta, sbIndexRead, r1, r2, sbGroupMap)
    }

    TAG_DNA_SAMPLE_BARCODE(ch_sb_input, tagHelperScripts, tagCodonScripts)
    ch_versions = ch_versions.mix(TAG_DNA_SAMPLE_BARCODE.out.versions)

    // Add modality barcodes from the tagmentation-specific DNA index stream.
    ch_mo_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        def moIndexRead = selectDnaIndexRead(meta, i1, i2, 'dna_modality_index_read')
        tuple(sampleId, meta, moIndexRead, modalityWhitelist)
    }

    ch_mo_input = ch_mo_meta
        .join(TAG_DNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, indexRead, modalityWhitelist, metaFromTag, taggedR1, taggedR2, readSetCounts ->
            tuple(sampleId, metaFromInput, indexRead, taggedR1, taggedR2, modalityWhitelist, readSetCounts)
        }

    TAG_DNA_MODALITY_BARCODE(ch_mo_input, tagHelperScripts, tagCodonScripts)
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
        .map { sampleId, metaFromInput, ligationRead, cellWhitelist, metaFromTag, taggedR1, taggedR2, readSetCounts ->
            tuple(sampleId, metaFromInput, ligationRead, taggedR1, taggedR2, cellWhitelist, readSetCounts)
        }

    TAG_DNA_CELL_BARCODE(ch_cb_input, cellHelperScripts, cellCodonScripts)
    ch_versions = ch_versions.mix(TAG_DNA_CELL_BARCODE.out.versions)

    // Every DNA pair first undergoes ordinary adapter/quality trimming so
    // genomic prefixes can be salvaged before exact residual-linker detection.
    TRIM_DNA_FASTQS(TAG_DNA_CELL_BARCODE.out.tagged, trimHelperScript)
    ch_versions = ch_versions.mix(TRIM_DNA_FASTQS.out.versions)

    // Only dual-tagmentation DNA is eligible for exact linker-signature
    // cleanup. A branch operation gives every trimmed pair to exactly one
    // route; single-tagmentation and explicitly disabled samples never invoke
    // the Cutadapt process. Both routes are recombined once before splitting.
    ch_dual_tag_artifact_routes = TRIM_DNA_FASTQS.out.trimmed.branch { _sampleId, meta, _trimmedR1, _trimmedR2 ->
        filter: params.filter_dual_tag_artifacts && meta.dna_tagmentation == 'dual'
        bypass: true
    }

    DUAL_TAG_ARTIFACT_FILTER(
        ch_dual_tag_artifact_routes.filter,
        file("${projectDir}/assets/dual_tag_artifact_23mers.fasta"),
        dualTagFilterHelperScript
    )
    ch_versions = ch_versions.mix(DUAL_TAG_ARTIFACT_FILTER.out.versions)

    ch_trimmed_for_split = ch_dual_tag_artifact_routes.bypass
        .mix(DUAL_TAG_ARTIFACT_FILTER.out.filtered)

    // Count the exact saved barcode decisions only for pairs that reach
    // splitting. For filtered dual-tag libraries this is the post-artifact
    // population; for every other DNA library it is the paired-trim output.
    ch_dna_tag_records_by_sample = TAG_DNA_CELL_BARCODE.out.tres_tag_records
        .map { tagRecords ->
            def sampleId = tagRecords.getName().replaceFirst(/\.dna_tag_records\.tsv\.gz$/, '')
            tuple(sampleId, tagRecords)
        }
    ch_dna_barcode_metric_maps = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, moMap, sbGroupMap)
    }
    ch_dna_barcode_metric_input = ch_trimmed_for_split
        .join(ch_dna_tag_records_by_sample)
        .join(ch_dna_barcode_metric_maps)
        .map { sampleId, splitMeta, splitR1, splitR2, tagRecords, contractMeta, moMap, sbGroupMap ->
            tuple(sampleId, splitMeta, 'dna', splitR1, splitR2, tagRecords, [sbGroupMap, moMap])
        }

    def barcodeMetricHelpers = [
        file("${projectDir}/bin/write_barcode_gate_metrics.py", checkIfExists: true),
        file("${projectDir}/bin/run_split_reads_dna.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    BARCODE_GATE_METRICS(ch_dna_barcode_metric_input, barcodeMetricHelpers)
    ch_versions = ch_versions.mix(BARCODE_GATE_METRICS.out.versions)

    // Split post-trimming DNA reads by sample-barcode group and modality mark.
    ch_split_meta = ch_dna_samples.map { sampleId, meta, i1, i2, r1, r2, modalityWhitelist, cellWhitelist, moMap, sbGroupMap ->
        tuple(sampleId, meta, moMap, sbGroupMap)
    }

    ch_split_input = ch_split_meta
        .join(ch_trimmed_for_split)
        .map { sampleId, metaFromInput, moMap, sbGroupMap, metaFromTrim, trimmedR1, trimmedR2 ->
            tuple(sampleId, metaFromInput, trimmedR1, trimmedR2, moMap, sbGroupMap)
        }

    SPLIT_DNA_READS(ch_split_input, splitHelperScripts, splitCodonScripts)
    ch_versions = ch_versions.mix(SPLIT_DNA_READS.out.versions)

    // The plain split FASTQs branch independently: alignment always consumes
    // them directly. The compression process runs only when split publication
    // is enabled.
    ch_published_split_fastqs = channel.empty()
    if( params.publish_split_fastqs ) {
        ch_compress_split_input = SPLIT_DNA_READS.out.split_fastqs
            .map { sampleId, meta, splitR1s, splitR2s ->
                tuple(sampleId, meta, 'dna', splitR1s, splitR2s)
            }

        COMPRESS_DNA_SPLIT_FASTQS(ch_compress_split_input)
        ch_versions = ch_versions.mix(COMPRESS_DNA_SPLIT_FASTQS.out.versions)
        ch_published_split_fastqs = COMPRESS_DNA_SPLIT_FASTQS.out.compressed_fastqs
    }

    ch_align_fastqs = SPLIT_DNA_READS.out.split_fastqs
        .flatMap { sampleId, meta, splitR1s, splitR2s ->
            pairDnaSplitFastqs(splitR1s, splitR2s).collect { split ->
                tuple(split.splitName, sampleId, meta, split.r1, split.r2)
            }
        }

    ch_align_rg = SPLIT_DNA_READS.out.rg_headers
        .flatMap { sampleId, meta, rgHeaders ->
            collectDnaRgHeaders(rgHeaders).collect { rg ->
                tuple(rg.splitName, sampleId, meta, rg.rgHeader)
            }
        }

    ch_align_input = ch_align_fastqs
        .join(ch_align_rg)
        .map { splitName, sampleId, metaFromFastq, splitR1, splitR2, sampleIdFromRg, metaFromRg, rgHeader ->
            def splitMeta = parseDnaSplitName(sampleId, splitName)
            def bwaReference = metaFromFastq.dna_bwa_reference as String
            def bwaReferenceName = new File(bwaReference).name
            def bwaIndexFiles = bwaMem2IndexSuffixes.collect { suffix ->
                file("${bwaReference}${suffix}", checkIfExists: true)
            }

            tuple(
                splitName,
                metaFromFastq,
                splitMeta.sampleGroup,
                splitMeta.modality,
                splitR1,
                splitR2,
                rgHeader,
                bwaReferenceName,
                bwaIndexFiles,
                file(metaFromFastq.dna_blacklist_bed as String, checkIfExists: true),
                (metaFromFastq.dna_effective_genome_size as String)
            )
        }

    // Finish the DNA core with alignment, duplicate marking, NoDup extraction, and coverage.
    ALIGN_DNA(ch_align_input, alignRuntimeScripts)
    ch_versions = ch_versions.mix(ALIGN_DNA.out.versions)

    // Preserve MarkDuplicates input exactly as before. A separate canonical
    // copy supplies the aligned-BAM output and QC branch.
    ch_canonical_aligned_input = ALIGN_DNA.out.bam
        .join(ALIGN_DNA.out.bai)
        .map { splitName, metaFromBam, alignedBam, metaFromBai, alignedBai ->
            tuple(
                splitName,
                metaFromBam,
                alignedBam,
                alignedBai,
                file(metaFromBam.canonical_chromosomes)
            )
        }

    FILTER_CANONICAL_DNA_ALIGNED_BAM(ch_canonical_aligned_input, canonicalBamRuntimeScripts)
    ch_versions = ch_versions.mix(FILTER_CANONICAL_DNA_ALIGNED_BAM.out.versions)

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
            def restoredMeta = restoreDnaMeta(nfMeta)
            tuple(
                nfMeta.tres_split_name as String,
                restoredMeta,
                markedDupBam,
                markedDupBai,
                markedDupMetrics,
                file(restoredMeta.canonical_chromosomes)
            )
        }

    NORMALIZE_DNA_MARKDUPLICATES(ch_normalize_markduplicates_input, canonicalBamRuntimeScripts)
    ch_versions = ch_versions.mix(NORMALIZE_DNA_MARKDUPLICATES.out.versions)

    ch_split_duplicates_input = NORMALIZE_DNA_MARKDUPLICATES.out.bam
        .map { splitName, meta, markedDupBam ->
            tuple(
                splitName,
                meta,
                markedDupBam,
                (meta.dna_effective_genome_size as String)
            )
        }

    SPLIT_DUPLICATES_DNA(ch_split_duplicates_input, splitDuplicatesRuntimeScripts)
    ch_versions = ch_versions.mix(SPLIT_DUPLICATES_DNA.out.versions)

    ch_nodup_for_coverage = SPLIT_DUPLICATES_DNA.out.bam
        .join(SPLIT_DUPLICATES_DNA.out.bai)
        .join(SPLIT_DUPLICATES_DNA.out.mapped_reads)
        .filter { splitName, metaFromBam, noDupBam, metaFromBai, noDupBai, metaFromCount, effectiveGenomeSize, mappedReadsFile ->
            hasMappedDnaReads(mappedReadsFile)
        }
        .map { splitName, metaFromBam, noDupBam, metaFromBai, noDupBai, metaFromCount, effectiveGenomeSize, mappedReadsFile ->
            def coverageMeta = metaFromBam + [
                dna_effective_genome_size: effectiveGenomeSize,
            ]
            tuple(
                nfcoreDnaMeta(splitName, coverageMeta, 'nodup_coverage'),
                noDupBam,
                noDupBai
            )
        }

    DEEPTOOLS_BAMCOVERAGE(
        ch_nodup_for_coverage,
        [],
        [],
        tuple([id: 'no_blacklist'], [])
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions_deeptools)
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions_samtools)

    ch_coverage_bigwigs = DEEPTOOLS_BAMCOVERAGE.out.bigwig.map { nfMeta, bigwig ->
        tuple(
            nfMeta.tres_split_name as String,
            restoreDnaMeta(nfMeta),
            bigwig
        )
    }

    ch_barcode_reports = TAG_DNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_DNA_MODALITY_BARCODE.out.metrics)
        .mix(TAG_DNA_CELL_BARCODE.out.metrics)

    ch_barcode_report_files = TAG_DNA_SAMPLE_BARCODE.out.metrics
        .flatMap { sampleId, counts, stats -> [counts, stats] }
        .mix(TAG_DNA_MODALITY_BARCODE.out.metrics.flatMap { sampleId, counts, stats -> [counts, stats] })
        .mix(TAG_DNA_CELL_BARCODE.out.metrics.flatMap { sampleId, counts, statsL1, statsL2, statsL3 ->
            [counts, statsL1, statsL2, statsL3]
        })
        .mix(BARCODE_GATE_METRICS.out.metrics.flatMap { sampleId, meta, gates, composition ->
            [gates, composition]
        })

    emit:
    tagged_fastqs   = TAG_DNA_CELL_BARCODE.out.tagged
    dual_tag_artifact_filter_qc = DUAL_TAG_ARTIFACT_FILTER.out.qc
    trimmed_fastqs  = TRIM_DNA_FASTQS.out.trimmed
    split_fastqs    = ch_published_split_fastqs
    rg_headers      = SPLIT_DNA_READS.out.rg_headers
    split_retention_metrics = SPLIT_DNA_READS.out.retention_metrics
    barcode_gate_metrics = BARCODE_GATE_METRICS.out.metrics
    alignment_retention_metrics = ALIGN_DNA.out.retention_metrics
    aligned_bams    = FILTER_CANONICAL_DNA_ALIGNED_BAM.out.bam
    aligned_bais    = FILTER_CANONICAL_DNA_ALIGNED_BAM.out.bai
    markeddup_bams = NORMALIZE_DNA_MARKDUPLICATES.out.bam
    markeddup_bais = NORMALIZE_DNA_MARKDUPLICATES.out.bai
    duplicate_metrics = NORMALIZE_DNA_MARKDUPLICATES.out.metrics
    nodup_bams = SPLIT_DUPLICATES_DNA.out.bam
    nodup_bais = SPLIT_DUPLICATES_DNA.out.bai
    coverage_bigwigs = ch_coverage_bigwigs
    coverage_warnings = SPLIT_DUPLICATES_DNA.out.warnings
    barcode_reports = ch_barcode_reports
    barcode_report_files = ch_barcode_report_files
    tres_tag_records = TAG_DNA_CELL_BARCODE.out.tres_tag_records
    versions = ch_versions
}
