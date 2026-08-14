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
 *      steps plus DNA trim_galore, dual-tag artifact filtering, Split_ReadsV2 dna mode, AlignDNA.sh,
 *      GATK MarkDuplicates, duplicate filtering to NoDup BAMs, and bamCoverage.
 */

include { RNA_CORE } from '../subworkflows/local/rna_core'
include { DNA_CORE } from '../subworkflows/local/dna_core'
include { TRES_REPORT_HTML } from '../modules/local/tres_report_html/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { SAMTOOLS_BAM_QC } from '../modules/local/samtools_bam_qc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'

def toRnaCoreInput(row) {
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

def toDnaCoreInput(row) {
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
        barcode_defaults          : params.barcode_defaults,
    ]
}

def qcMeta(meta, id, modality, stage, splitName) {
    return meta + [
        id              : id,
        tres_modality   : modality,
        tres_qc_stage   : stage,
        tres_split_name : splitName,
    ]
}

def toFastqcInput(row) {
    def reads = row.modality == 'dna'
        ? [row.i1, row.i2, row.r1, row.r2]
        : [row.i1, row.r1, row.r2]

    return tuple(
        qcMeta(row, "${row.modality}.${row.id}.raw", row.modality as String, 'raw_fastq', row.id as String),
        reads.findAll { it }.collect { file(it) }
    )
}

def validateCoreResourceContract(rnaRows, dnaRows, maxCpus) {
    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }
}

workflow TRESEQ {
    take:
    sampleRows

    main:
    // Parse the single supported YAML contract into modality-specific work rows.
    def rnaRows = sampleRows.findAll { row -> row.modality == 'rna' }
    def dnaRows = sampleRows.findAll { row -> row.modality == 'dna' }
    def maxCpus = params.max_cpus as int
    def libraryName = sampleRows ? (sampleRows[0].library_name as String) : 'unknown library'

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

    Channel
        .fromList(sampleRows)
        .map { row -> toFastqcInput(row) }
        .set { ch_raw_fastqs_for_fastqc }

    FASTQC(ch_raw_fastqs_for_fastqc)

    // Combined Samtools sidecar QC does not alter the TrESFlow data path; one
    // task per BAM emits the standardized QC text files.
    ch_rna_bams_for_qc = RNA_CORE.out.internal_filtered_bams.map { splitName, meta, bam ->
        tuple(qcMeta(meta, "rna.${splitName}.filtered_cells", 'rna', 'filtered_cells', splitName), bam, [], false)
    }

    ch_dna_aligned_bams_for_qc = DNA_CORE.out.aligned_bams
        .join(DNA_CORE.out.aligned_bais)
        .map { splitName, metaFromBam, bam, metaFromBai, bai ->
            tuple(qcMeta(metaFromBam, "dna.${splitName}.aligned", 'dna', 'aligned', splitName), bam, bai, true)
        }

    ch_dna_markeddup_bams_for_qc = DNA_CORE.out.markeddup_bams
        .join(DNA_CORE.out.markeddup_bais)
        .map { splitName, metaFromBam, bam, metaFromBai, bai ->
            tuple(qcMeta(metaFromBam, "dna.${splitName}.markeddup", 'dna', 'markeddup', splitName), bam, bai, true)
        }

    ch_dna_nodup_bams_for_qc = DNA_CORE.out.nodup_bams
        .join(DNA_CORE.out.nodup_bais)
        .map { splitName, metaFromBam, bam, metaFromBai, bai ->
            tuple(qcMeta(metaFromBam, "dna.${splitName}.nodup", 'dna', 'nodup', splitName), bam, bai, true)
        }

    ch_bams_for_samtools_qc = ch_rna_bams_for_qc
        .mix(ch_dna_aligned_bams_for_qc)
        .mix(ch_dna_markeddup_bams_for_qc)
        .mix(ch_dna_nodup_bams_for_qc)

    SAMTOOLS_BAM_QC(ch_bams_for_samtools_qc)

    ch_barcode_report_files = RNA_CORE.out.barcode_report_files
        .mix(DNA_CORE.out.barcode_report_files)

    ch_report_source_files = ch_barcode_report_files
        .mix(DNA_CORE.out.dual_tag_artifact_filter_qc.flatMap { _sampleId, _meta, cutadaptJson, summary ->
            [cutadaptJson, summary]
        })
        .mix(RNA_CORE.out.aligned_solo_summaries.map { splitName, meta, soloSummary -> soloSummary })
        .mix(RNA_CORE.out.aligned_star_logs.map { splitName, meta, starLog -> starLog })
        .mix(DNA_CORE.out.duplicate_metrics.map { splitName, meta, metrics -> metrics })
        .mix(FASTQC.out.zip.map { meta, zip -> zip })
        .mix(SAMTOOLS_BAM_QC.out.flagstat.map { meta, flagstat -> flagstat })
        .mix(SAMTOOLS_BAM_QC.out.stats.map { meta, stats -> stats })
        .mix(SAMTOOLS_BAM_QC.out.idxstats.map { meta, idxstats -> idxstats })
        .mix(SAMTOOLS_BAM_QC.out.quickcheck.map { meta, report -> report })

    ch_tres_report_input = ch_report_source_files
        .collect()
        .map { files -> tuple([id: 'tresflow', library_name: libraryName], files) }

    ch_multiqc_input = ch_report_source_files
        .collect()
        .map { files ->
            tuple(
                [id: 'tresflow'],
                files,
                file("${projectDir}/assets/multiqc_config.yml"),
                [],
                [],
                []
            )
        }

    TRES_REPORT_HTML(ch_tres_report_input)
    MULTIQC(ch_multiqc_input)

    emit:
    tagged_fastqs               = RNA_CORE.out.tagged_fastqs
    trimmed_fastqs              = RNA_CORE.out.trimmed_fastqs
    split_fastqs                = RNA_CORE.out.split_fastqs
    rg_headers                  = RNA_CORE.out.rg_headers
    rna_split_retention_metrics = RNA_CORE.out.split_retention_metrics
    rna_barcode_gate_metrics = RNA_CORE.out.barcode_gate_metrics
    rna_filter_retention_metrics = RNA_CORE.out.filter_retention_metrics
    usam_files                  = RNA_CORE.out.usam_files
    aligned_solo_dirs           = RNA_CORE.out.aligned_solo_dirs
    aligned_solo_summaries      = RNA_CORE.out.aligned_solo_summaries
    aligned_star_logs           = RNA_CORE.out.aligned_star_logs
    aligned_filtered_bams       = RNA_CORE.out.aligned_filtered_bams
    aligned_stranded_bigwigs    = RNA_CORE.out.aligned_stranded_bigwigs
    aligned_unstranded_bigwigs = RNA_CORE.out.aligned_unstranded_bigwigs
    barcode_reports             = RNA_CORE.out.barcode_reports
    dna_tagged_fastqs           = DNA_CORE.out.tagged_fastqs
    dna_dual_tag_artifact_filter_qc = DNA_CORE.out.dual_tag_artifact_filter_qc
    dna_trimmed_fastqs          = DNA_CORE.out.trimmed_fastqs
    dna_split_fastqs            = DNA_CORE.out.split_fastqs
    dna_rg_headers              = DNA_CORE.out.rg_headers
    dna_split_retention_metrics = DNA_CORE.out.split_retention_metrics
    dna_barcode_gate_metrics = DNA_CORE.out.barcode_gate_metrics
    dna_alignment_retention_metrics = DNA_CORE.out.alignment_retention_metrics
    dna_aligned_bams            = DNA_CORE.out.aligned_bams
    dna_aligned_bais            = DNA_CORE.out.aligned_bais
    dna_markeddup_bams          = DNA_CORE.out.markeddup_bams
    dna_markeddup_bais          = DNA_CORE.out.markeddup_bais
    dna_duplicate_metrics       = DNA_CORE.out.duplicate_metrics
    dna_nodup_bams              = DNA_CORE.out.nodup_bams
    dna_nodup_bais              = DNA_CORE.out.nodup_bais
    dna_coverage_bigwigs        = DNA_CORE.out.coverage_bigwigs
    dna_coverage_warnings       = DNA_CORE.out.coverage_warnings
    dna_barcode_reports         = DNA_CORE.out.barcode_reports
    samtools_flagstat           = SAMTOOLS_BAM_QC.out.flagstat
    samtools_stats              = SAMTOOLS_BAM_QC.out.stats
    samtools_idxstats           = SAMTOOLS_BAM_QC.out.idxstats
    samtools_quickcheck         = SAMTOOLS_BAM_QC.out.quickcheck
    fastqc_html                 = FASTQC.out.html
    fastqc_zip                  = FASTQC.out.zip
    multiqc_report              = MULTIQC.out.report
    multiqc_data                = MULTIQC.out.data
    tres_report_html            = TRES_REPORT_HTML.out.html
    tres_report_metrics_json    = TRES_REPORT_HTML.out.metrics_json
}
