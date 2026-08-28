/*
 * Subworkflow: RNA_CORE
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - ordered raw RNA I1 / R1 / R2 FASTQ collections (one logical library row)
 *   - RNA cell-barcode whitelist
 *   - internally derived sample-barcode group map TSV, used both to derive the
 *     effective RNA sample-barcode whitelist and to split grouped RNA reads
 *   - exact RNA STAR index directory from references.rna_ref_dir and derived chromosome sizes carried through sample metadata
 * Outputs:
 *   - RNA FASTQs tagged with SB, UM, then CB comments
 *   - uncompressed trim_galore paired-end FASTQs from the CB-tagged reads
 *   - plain Split_ReadsV2 per-group RNA FASTQs for computation plus optional independently compressed publication copies and SAM RG headers
 *   - FqToSAM unmapped SAM files from each split RNA FASTQ pair
 *   - STARsolo outputs, low-compression published filtered BAMs, and bigWigs
 *   - barcode count/stat files from all wrapped RNA steps
 */

include { TAG_RNA_SAMPLE_BARCODE } from '../../../modules/local/tag_rna_sb/main'
include { TAG_RNA_UMI }            from '../../../modules/local/tag_rna_umi/main'
include { TAG_RNA_CELL_BARCODE }   from '../../../modules/local/tag_rna_cell_barcode/main'
include { TRIM_RNA_FASTQS }        from '../../../modules/local/trim_rna_fastqs/main'
include { BARCODE_GATE_METRICS }    from '../../../modules/local/barcode_gate_metrics/main'
include { SPLIT_RNA_READS }        from '../../../modules/local/split_rna_reads/main'
include { COMPRESS_SPLIT_FASTQS as COMPRESS_RNA_SPLIT_FASTQS } from '../../../modules/local/compress_split_fastqs/main'
include { FQ_TO_SAM }              from '../../../modules/local/fq_to_sam/main'
include { RNA_STARSOLO_ALIGN }     from '../../../modules/local/rna_starsolo_align/main'
include { RNA_FILTERED_BAM }       from '../../../modules/local/rna_filtered_bam/main'
include { RNA_COVERAGE }           from '../../../modules/local/rna_coverage/main'

def asRnaPathList(value) {
    return value instanceof List ? value : [value]
}

def pairRnaSplitFastqs(sampleId, splitR1s, splitR2s) {
    def r1ByGroup = asRnaPathList(splitR1s).collectEntries { path ->
        def group = path.getName()
            .replaceFirst("^${sampleId}_", '')
            .replaceFirst('_R1\\.(?:fastq|fq)(?:\\.gz)?$', '')
        [(group): path]
    }
    def r2ByGroup = asRnaPathList(splitR2s).collectEntries { path ->
        def group = path.getName()
            .replaceFirst("^${sampleId}_", '')
            .replaceFirst('_R2\\.(?:fastq|fq)(?:\\.gz)?$', '')
        [(group): path]
    }

    def groups = (r1ByGroup.keySet() + r2ByGroup.keySet()).unique().sort()
    groups.collect { group ->
        if( !r1ByGroup.containsKey(group) || !r2ByGroup.containsKey(group) ) {
            throw new IllegalStateException(
                "Missing split FASTQ mate for sample '${sampleId}' group '${group}'"
            )
        }
        [splitName: "${sampleId}_${group}", r1: r1ByGroup[group], r2: r2ByGroup[group]]
    }
}

workflow RNA_CORE {
    take:
    ch_rna_samples

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
    def umiHelperScripts = [
        file("${projectDir}/bin/run_tag_umi.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def umiCodonScripts = [
        file("${codonScriptsRoot}/Tag_UMI.codon", checkIfExists: true),
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
        file("${projectDir}/bin/run_split_reads_rna.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def splitCodonScripts = [
        file("${codonScriptsRoot}/Split_ReadsV2.codon", checkIfExists: true),
    ]
    def fqToSamHelperScripts = [
        file("${projectDir}/bin/run_fq_to_sam.py", checkIfExists: true),
    ]
    def fqToSamCodonScripts = [
        file("${codonScriptsRoot}/FqToSAM.codon", checkIfExists: true),
    ]

    // Tag sample barcodes from the RNA read-2 stream.
    ch_sb_input = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, r1, r2, sbGroupMap)
    }

    TAG_RNA_SAMPLE_BARCODE(ch_sb_input, tagHelperScripts, tagCodonScripts)
    ch_versions = ch_versions.mix(TAG_RNA_SAMPLE_BARCODE.out.versions)

    // Add UMIs after sample-barcode tagging.
    ch_raw_r2 = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, r2)
    }

    ch_umi_input = ch_raw_r2
        .join(TAG_RNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, rawR2, metaFromTag, taggedR1, taggedR2, readSetCounts ->
            tuple(sampleId, metaFromInput, rawR2, taggedR1, taggedR2, readSetCounts)
        }

    TAG_RNA_UMI(ch_umi_input, umiHelperScripts, umiCodonScripts)
    ch_versions = ch_versions.mix(TAG_RNA_UMI.out.versions)

    // Add ligation-derived cell barcodes from I1.
    ch_cb_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, i1, cellWhitelist)
    }

    ch_cb_input = ch_cb_meta
        .join(TAG_RNA_UMI.out.tagged)
        .map { sampleId, metaFromInput, i1, cellWhitelist, metaFromTag, taggedR1, taggedR2, readSetCounts ->
            tuple(sampleId, metaFromInput, i1, taggedR1, taggedR2, cellWhitelist, readSetCounts)
        }

    TAG_RNA_CELL_BARCODE(ch_cb_input, cellHelperScripts, cellCodonScripts)
    ch_versions = ch_versions.mix(TAG_RNA_CELL_BARCODE.out.versions)
    TRIM_RNA_FASTQS(TAG_RNA_CELL_BARCODE.out.tagged)
    ch_versions = ch_versions.mix(TRIM_RNA_FASTQS.out.versions)

    // Count exact same-pair cumulative gates from the barcode decisions
    // already saved by Tag_Lig3, restricted to pairs surviving paired
    // trimming. This sidecar has no FASTQ output and cannot change routing.
    ch_rna_tag_records_by_sample = TAG_RNA_CELL_BARCODE.out.tres_tag_records
        .map { tagRecords ->
            def sampleId = tagRecords.getName().replaceFirst(/\.rna_tag_records\.tsv\.gz$/, '')
            tuple(sampleId, tagRecords)
        }
    ch_rna_barcode_metric_maps = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, sbGroupMap)
    }
    ch_rna_barcode_metric_input = TRIM_RNA_FASTQS.out.trimmed
        .join(ch_rna_tag_records_by_sample)
        .join(ch_rna_barcode_metric_maps)
        .map { sampleId, trimMeta, trimmedR1, trimmedR2, tagRecords, contractMeta, sbGroupMap ->
            tuple(sampleId, trimMeta, 'rna', trimmedR1, trimmedR2, tagRecords, [sbGroupMap])
        }

    def barcodeMetricHelpers = [
        file("${projectDir}/bin/write_barcode_gate_metrics.py", checkIfExists: true),
        file("${projectDir}/bin/run_split_reads_dna.py", checkIfExists: true),
        file("${projectDir}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    BARCODE_GATE_METRICS(ch_rna_barcode_metric_input, barcodeMetricHelpers)
    ch_versions = ch_versions.mix(BARCODE_GATE_METRICS.out.versions)

    // Split trimmed reads by sample-barcode group before FQ_TO_SAM.
    ch_split_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, sbGroupMap)
    }

    ch_split_input = ch_split_meta
        .join(TRIM_RNA_FASTQS.out.trimmed)
        .map { sampleId, metaFromInput, sbGroupMap, metaFromTrim, trimmedR1, trimmedR2 ->
            tuple(sampleId, metaFromInput, trimmedR1, trimmedR2, sbGroupMap)
        }

    SPLIT_RNA_READS(ch_split_input, splitHelperScripts, splitCodonScripts)
    ch_versions = ch_versions.mix(SPLIT_RNA_READS.out.versions)

    // The plain split FASTQs branch independently: computation always consumes
    // them directly. The compression process runs only when split publication
    // is enabled.
    ch_published_split_fastqs = channel.empty()
    if( params.publish_split_fastqs ) {
        ch_compress_split_input = SPLIT_RNA_READS.out.split_fastqs
            .map { sampleId, meta, splitR1s, splitR2s ->
                tuple(sampleId, meta, 'rna', splitR1s, splitR2s)
            }

        COMPRESS_RNA_SPLIT_FASTQS(ch_compress_split_input)
        ch_versions = ch_versions.mix(COMPRESS_RNA_SPLIT_FASTQS.out.versions)
        ch_published_split_fastqs = COMPRESS_RNA_SPLIT_FASTQS.out.compressed_fastqs
    }

    ch_fq_to_sam_input = SPLIT_RNA_READS.out.split_fastqs
        .flatMap { sampleId, meta, splitR1s, splitR2s ->
            pairRnaSplitFastqs(sampleId, splitR1s, splitR2s).collect { split ->
                tuple(split.splitName, meta, split.r1, split.r2)
            }
        }

    FQ_TO_SAM(ch_fq_to_sam_input, fqToSamHelperScripts, fqToSamCodonScripts)
    ch_versions = ch_versions.mix(FQ_TO_SAM.out.versions)

    // Align each grouped unmapped SAM independently with STARsolo.
    ch_starsolo_input = FQ_TO_SAM.out.usam
        .map { splitName, meta, usam ->
            tuple(
                splitName,
                meta,
                usam,
                meta.rna_star_index_dir as String
            )
        }

    RNA_STARSOLO_ALIGN(ch_starsolo_input)
    ch_versions = ch_versions.mix(RNA_STARSOLO_ALIGN.out.versions)

    ch_filtered_bam_input = RNA_STARSOLO_ALIGN.out.solo_dir
        .join(RNA_STARSOLO_ALIGN.out.aligned_bam)
        .map { splitName, metaFromSolo, soloDir, metaFromBam, alignedBam ->
            tuple(
                splitName,
                metaFromSolo,
                soloDir,
                alignedBam,
                file(metaFromSolo.canonical_chromosomes)
            )
        }

    RNA_FILTERED_BAM(ch_filtered_bam_input)
    ch_versions = ch_versions.mix(RNA_FILTERED_BAM.out.versions)

    ch_coverage_input = RNA_FILTERED_BAM.out.filtered_bam
        .map { splitName, meta, filteredBam ->
            tuple(
                splitName,
                meta,
                filteredBam,
                meta.rna_star_index_dir as String,
                meta.canonical_chrom_sizes as String
            )
        }

    RNA_COVERAGE(ch_coverage_input)
    ch_versions = ch_versions.mix(RNA_COVERAGE.out.versions)

    ch_barcode_reports = TAG_RNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_RNA_UMI.out.metrics)
        .mix(TAG_RNA_CELL_BARCODE.out.metrics)

    // Use the modality-qualified publication copies for report consumers. The
    // generic runtime filenames are intentionally retained only on the
    // computational tuples and are ambiguous when RNA and DNA share a sample.
    ch_barcode_report_files = TAG_RNA_SAMPLE_BARCODE.out.tres_stats
        .flatMap { reportFiles -> asRnaPathList(reportFiles) }
        .mix(TAG_RNA_UMI.out.tres_stats.flatMap { reportFiles -> asRnaPathList(reportFiles) })
        .mix(TAG_RNA_CELL_BARCODE.out.tres_cell_stats.flatMap { reportFiles -> asRnaPathList(reportFiles) })
        .mix(BARCODE_GATE_METRICS.out.metrics.flatMap { sampleId, meta, gates, composition ->
            [gates, composition]
        })

    emit:
    tagged_fastqs    = TAG_RNA_CELL_BARCODE.out.tagged
    trimmed_fastqs   = TRIM_RNA_FASTQS.out.trimmed
    split_fastqs     = ch_published_split_fastqs
    rg_headers       = SPLIT_RNA_READS.out.rg_headers
    split_retention_metrics = SPLIT_RNA_READS.out.retention_metrics
    barcode_gate_metrics = BARCODE_GATE_METRICS.out.metrics
    usam_files       = FQ_TO_SAM.out.usam
    aligned_solo_dirs = RNA_STARSOLO_ALIGN.out.solo_dir
    aligned_solo_summaries = RNA_STARSOLO_ALIGN.out.solo_summary
    aligned_solo_report_summaries = RNA_STARSOLO_ALIGN.out.report_solo_summary
    aligned_star_logs = RNA_STARSOLO_ALIGN.out.star_log
    internal_filtered_bams = RNA_FILTERED_BAM.out.filtered_bam
    filter_retention_metrics = RNA_FILTERED_BAM.out.retention_metrics
    aligned_filtered_bams = RNA_FILTERED_BAM.out.filtered_bam
    aligned_stranded_bigwigs = RNA_COVERAGE.out.stranded_bw
    aligned_unstranded_bigwigs = RNA_COVERAGE.out.unstranded_bw
    barcode_reports  = ch_barcode_reports
    barcode_report_files = ch_barcode_report_files
    tres_tag_records = TAG_RNA_CELL_BARCODE.out.tres_tag_records
    versions         = ch_versions
}
