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

include { BULK_REPORT } from '../modules/local/bulk_report/main'
include { RNA_CORE } from '../subworkflows/local/rna_core'
include { DNA_CORE } from '../subworkflows/local/dna_core'
include { TRES_REPORT_HTML } from '../modules/local/tres_report_html/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { SAMTOOLS_BAM_QC } from '../modules/local/samtools_bam_qc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { VALIDATE_PREPROCESSED_SPLIT } from '../modules/local/validate_preprocessed_split/main'
include { WRITE_PREPROCESSED_MANIFEST } from '../modules/local/write_preprocessed_manifest/main'

def toRnaCoreInput(row) {
    tuple(
        row.id,
        row,
        row.i1.collect { path -> file(path) },
        row.r1.collect { path -> file(path) },
        row.r2.collect { path -> file(path) },
        file(row.cell_whitelist),
        file(row.sb_group_map)
    )
}

def toDnaCoreInput(row) {
    tuple(
        row.id,
        row,
        row.i1.collect { path -> file(path) },
        row.i2.collect { path -> file(path) },
        row.r1.collect { path -> file(path) },
        row.r2.collect { path -> file(path) },
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
    def readRoles = row.modality == 'dna'
        ? (row.i2_implicit ? [row.i1, row.r1, row.r2] : [row.i1, row.i2, row.r1, row.r2])
        : [row.i1, row.r1, row.r2]
    def reads = readRoles.collectMany { roleReads -> roleReads as List }

    return tuple(
        qcMeta(row, "${row.modality}.${row.id}.raw", row.modality as String, 'raw_fastq', row.id as String),
        reads.collect { file(it) }
    )
}

def validateCoreResourceContract(rnaRows, dnaRows, maxCpus) {
    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }
}

def manifestEntries(meta, reads) {
    def dna = meta.modality == 'dna'
    def directory = dna ? 'dna_split_fastqs' : 'restart/rna'
    (reads instanceof List ? reads : [reads]).collect { r1 ->
        def splitId = r1.name.replaceFirst(/_R1\.fastq\.gz$/, '')
        def target = meta.split_targets[splitId]
        def item = [
            split_id: splitId, sample: meta.id, library: meta.library_name,
            group: target.group, modality: meta.modality, format: 'tresflow-split-v1',
            sb_mapping_sha256: meta.sb_mapping_sha256,
            r1: "${directory}/${splitId}_R1.fastq.gz",
            r2: "${directory}/${splitId}_R2.fastq.gz",
        ]
        if( dna ) item += [
            mark: target.mark, dna_tagmentation: meta.dna_tagmentation,
            rg: "restart/dna/SAM_RG_Header_${splitId}.tsv",
        ]
        item
    }
}

workflow TRESEQ {
    take:
    sampleRows
    reportMetadata

    main:
    // Parse the single supported YAML contract into modality-specific work rows.
    def restart = reportMetadata.restart
    def rawRows = restart ? [] : sampleRows
    def rnaRows = rawRows.findAll { row -> row.modality == 'rna' }
    def dnaRows = rawRows.findAll { row -> row.modality == 'dna' }
    def maxCpus = params.max_cpus as int

    validateCoreResourceContract(rnaRows, dnaRows, maxCpus)

    Channel
        .fromList(rnaRows)
        .map { row -> toRnaCoreInput(row) }
        .set { ch_rna_samples }

    Channel
        .fromList(dnaRows)
        .map { row -> toDnaCoreInput(row) }
        .set { ch_dna_samples }

    ch_restart_rna_pairs = Channel.empty()
    ch_restart_dna_pairs = Channel.empty()
    ch_routed_qc = Channel.empty()
    if( restart ) {
        ch_restart_sources = Channel.fromList(restart.splits)
        ch_validation = ch_restart_sources.map { split ->
            def entry = split.entry
            def meta = split.row
            def spec = entry + [
                bulk: reportMetadata.bulk,
                fastq_only: restart.fastq_only == true,
                r1: split.r1.toString(),
                r2: split.r2.toString(),
                sb_indices: (meta.cell_identity_records ?: []).findAll { it.group == entry.group }
                    .collectEntries { [(it.sb_bc): it.sb_index] },
                ligation_length: 3 * meta.cell_bc_len,
                umi_start: meta.umi_bc_start,
                umi_length: meta.umi_bc_len,
                mo_barcode: meta.split_targets[entry.split_id].mo_barcode,
            ]
            tuple(entry.split_id, meta, spec, split.r1, split.r2, split.rg ?: [])
        }
        // Bulk RNA validates and converts in the same pass over the supplied pair.
        def bulk = reportMetadata.bulk
        VALIDATE_PREPROCESSED_SPLIT(ch_validation.filter { id, meta, spec, r1, r2, rg -> !(bulk && meta.modality == 'rna') })
        if( bulk ) {
            ch_routed_qc = VALIDATE_PREPROCESSED_SPLIT.out.validated.map { id, meta, validated ->
                validated.resolve("${id}.dna.routed.json")
            }
        }
        ch_validated_pairs = VALIDATE_PREPROCESSED_SPLIT.out.validated
            .map { splitId, meta, validated -> tuple([meta.modality, splitId], meta, validated) }
            .join(ch_restart_sources.map { split ->
                tuple([split.entry.modality, split.entry.split_id], split.r1, split.r2)
            })
            .branch { key, meta, validated, r1, r2 ->
                rna: key[0] == 'rna'
                    tuple(key[1], meta + [group: meta.split_targets[key[1]].group, restart_input: true], r1, r2)
                dna: key[0] == 'dna'
                    tuple(key[1], meta + [restart_input: true], r1, r2, validated.resolve("SAM_RG_Header_${key[1]}.tsv"))
            }
        ch_restart_rna_pairs = bulk ? ch_validation.filter { id, meta, spec, r1, r2, rg -> meta.modality == 'rna' }
            .map { id, meta, spec, r1, r2, rg ->
                tuple(id, meta + [group: spec.group, restart_input: true, bulk_validation: spec], r1, r2)
            } : ch_validated_pairs.rna
        ch_restart_dna_pairs = ch_validated_pairs.dna
    }

    RNA_CORE(ch_rna_samples, reportMetadata.preprocess_only, ch_restart_rna_pairs, !!restart, reportMetadata.bulk)
    DNA_CORE(ch_dna_samples, reportMetadata.preprocess_only, ch_restart_dna_pairs, !!restart, reportMetadata.bulk)

    Channel
        .fromList(rawRows)
        .map { row -> toFastqcInput(row) }
        .set { ch_raw_fastqs_for_fastqc }

    FASTQC(ch_raw_fastqs_for_fastqc)

    ch_samtools_flagstat = Channel.empty()
    ch_samtools_stats = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    ch_samtools_quickcheck = Channel.empty()
    if( !reportMetadata.preprocess_only && !reportMetadata.bulk ) {
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

        ch_samtools_flagstat = SAMTOOLS_BAM_QC.out.flagstat
        ch_samtools_stats = SAMTOOLS_BAM_QC.out.stats
        ch_samtools_idxstats = SAMTOOLS_BAM_QC.out.idxstats
        ch_samtools_quickcheck = SAMTOOLS_BAM_QC.out.quickcheck
    }

    ch_barcode_report_files = RNA_CORE.out.barcode_report_files
        .mix(DNA_CORE.out.barcode_report_files)

    // MultiQC retains its broad source collection. The TrESFlow report instead
    // receives explicit metric contracts and waits on every producer, avoiding
    // output-directory scans and publication races.
    ch_multiqc_preprocessing_files = ch_barcode_report_files
        .mix(DNA_CORE.out.dual_tag_artifact_filter_qc.flatMap { _sampleId, _meta, cutadaptJson, summary ->
            [cutadaptJson, summary]
        })
        .mix(FASTQC.out.zip.flatMap { _meta, zips -> [zips].flatten() })

    ch_multiqc_source_files = ch_multiqc_preprocessing_files
        .mix(Channel.fromList(restart?.multiqc_inputs ?: []))
        .mix(RNA_CORE.out.aligned_solo_summaries.map { splitName, meta, soloSummary -> soloSummary })
        .mix(RNA_CORE.out.aligned_star_logs.map { splitName, meta, starLog -> starLog })
        .mix(DNA_CORE.out.duplicate_metrics.map { splitName, meta, metrics -> metrics })
        .mix(ch_samtools_flagstat.map { meta, flagstat -> flagstat })
        .mix(ch_samtools_stats.map { meta, stats -> stats })
        .mix(ch_samtools_idxstats.map { meta, idxstats -> idxstats })
        .mix(ch_samtools_quickcheck.map { meta, report -> report })

    ch_report_barcode_files = ch_barcode_report_files.filter { reportFile ->
        def name = reportFile.getName()
        name.endsWith('.stats.tsv') || name.endsWith('_sample_barcode.counts.tsv') ||
            name.endsWith('_barcode_gates.tsv') || name.endsWith('_barcode_composition.tsv') ||
            name ==~ /.+\.dna_cell\.stats_L[123]\.tsv/
    }

    def reportContractPaths = sampleRows
        .collectMany { row -> [row.sb_group_map, row.mo_map].findAll { it } }
        .collect { path -> file(path) }
        .unique { path -> path.toString() }

    ch_report_contract_files = Channel.fromList(reportContractPaths)

    ch_report_preprocessing_files = ch_report_barcode_files
        .mix(RNA_CORE.out.split_retention_metrics.map { sampleId, meta, metrics -> metrics })
        .mix(RNA_CORE.out.export_retention_metrics.map { sampleId, meta, metrics -> metrics })
        .mix(DNA_CORE.out.split_retention_metrics.map { sampleId, meta, metrics -> metrics })
        .mix(DNA_CORE.out.dual_tag_artifact_filter_qc.map { _sampleId, _meta, _cutadaptJson, summary -> summary })

    if( reportMetadata.preprocess_only ) {
        ch_manifest_splits = RNA_CORE.out.restart_fastqs.mix(DNA_CORE.out.split_fastqs)
            .flatMap { _sampleId, meta, r1, _r2 -> manifestEntries(meta, r1) }
            .toList()
        ch_manifest_headers = DNA_CORE.out.rg_headers
            .flatMap { _sampleId, _meta, headers -> [headers].flatten() }
            .toList()
        ch_manifest_metrics = ch_report_preprocessing_files.toList().map { files -> tuple('run', files) }
            .join(ch_multiqc_preprocessing_files.toList().map { files -> tuple('run', files) })
            .map { _run, reportFiles, multiqcFiles ->
                def files = (reportFiles + multiqcFiles).unique().sort { it.toString() }
                tuple(files, reportFiles.unique().collect { files.indexOf(it) }, multiqcFiles.unique().collect { files.indexOf(it) })
            }
        WRITE_PREPROCESSED_MANIFEST(
            [pipeline_version: reportMetadata.pipeline_version,
             runtime_env_prefix: sampleRows[0].runtime_env_prefix, runtime_tmpdir: sampleRows[0].runtime_tmpdir],
            ch_manifest_splits,
            ch_manifest_headers,
            ch_manifest_metrics
        )
    }

    ch_report_source_files = ch_report_preprocessing_files
        .mix(Channel.fromList(restart?.report_inputs ?: []))
        .mix(RNA_CORE.out.filter_retention_metrics.map { splitName, meta, metrics -> metrics })
        .mix(DNA_CORE.out.alignment_retention_metrics.map { splitName, meta, metrics -> metrics })
        .mix(RNA_CORE.out.aligned_solo_report_summaries.map { splitName, meta, reportSoloSummary -> reportSoloSummary })
        .mix(RNA_CORE.out.aligned_star_logs.map { splitName, meta, starLog -> starLog })
        .mix(DNA_CORE.out.duplicate_metrics.map { splitName, meta, metrics -> metrics })
        .mix(ch_samtools_flagstat.map { meta, flagstat -> flagstat })
        .mix(ch_report_contract_files)

    ch_tres_report_input = ch_report_source_files
        .collect()
        .map { files -> tuple([
            id              : 'tresflow',
            preprocess_only : reportMetadata.preprocess_only,
            report_title    : reportMetadata.report_title,
            pipeline_version: reportMetadata.pipeline_version,
            preprocessed_splits: restart && !restart.report_inputs ? restart.splits.collect { split ->
                def entry = split.entry
                [sample: entry.sample, group: entry.group, modality: entry.modality,
                 mark: split.row.split_targets[entry.split_id].mark, split_id: entry.split_id]
            } : [],
            filter_dual_tag_artifacts: restart ? restart.filter_dual_tag_artifacts : params.filter_dual_tag_artifacts,
            runtime_env_prefix: sampleRows ? sampleRows[0].runtime_env_prefix : '',
            runtime_tmpdir  : sampleRows ? sampleRows[0].runtime_tmpdir : '',
            samples         : sampleRows.collect { row -> [
                id              : row.id,
                modality        : row.modality,
                dna_tagmentation: row.dna_tagmentation,
                groups          : row.samplesheet_groups ?: [],
            ] },
        ], files) }

    ch_multiqc_input = ch_multiqc_source_files
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

    if( reportMetadata.bulk ) {
        def targets = sampleRows.collectMany { row ->
            row.split_targets.collect { id, target ->
                [id: id, sample: row.id, group: target.group, mark: target.mark,
                 library: row.library_name, modality: row.modality]
            }
        }
        ch_bulk_files = RNA_CORE.out.bulk_counts.mix(RNA_CORE.out.bulk_qc, DNA_CORE.out.bulk_qc, ch_routed_qc)
            .mix(ch_barcode_report_files.filter { f -> f.name.endsWith('_barcode_composition.tsv') ||
                f.name.endsWith('_barcode_gates.tsv') || f.name.endsWith('_sample_barcode.stats.tsv') })
        BULK_REPORT(ch_bulk_files.toList().map { files ->
            tuple(sampleRows[0], [targets: targets, restart: !!restart], files)
        })
    } else {
        TRES_REPORT_HTML(ch_tres_report_input)
        MULTIQC(ch_multiqc_input)
    }

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
    samtools_flagstat           = ch_samtools_flagstat
    samtools_stats              = ch_samtools_stats
    samtools_idxstats           = ch_samtools_idxstats
    samtools_quickcheck         = ch_samtools_quickcheck
    fastqc_html                 = FASTQC.out.html
    fastqc_zip                  = FASTQC.out.zip
    multiqc_report              = reportMetadata.bulk ? Channel.empty() : MULTIQC.out.report
    multiqc_data                = reportMetadata.bulk ? Channel.empty() : MULTIQC.out.data
    tres_report_html            = reportMetadata.bulk ? Channel.empty() : TRES_REPORT_HTML.out.html
    tres_report_read_retention  = reportMetadata.bulk ? Channel.empty() : TRES_REPORT_HTML.out.read_retention
    tres_report_qc_metrics      = reportMetadata.bulk ? Channel.empty() : TRES_REPORT_HTML.out.qc_metrics
    tres_report_barcode_composition = reportMetadata.bulk ? Channel.empty() : TRES_REPORT_HTML.out.barcode_composition
    tres_report_library_complexity = reportMetadata.bulk ? Channel.empty() : TRES_REPORT_HTML.out.library_complexity
}
