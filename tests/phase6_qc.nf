#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC } from '../modules/nf-core/fastqc/main'
include { SAMTOOLS_BAM_QC } from '../modules/local/samtools_bam_qc/main'
include { TRES_REPORT_HTML } from '../modules/local/tres_report_html/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'

def requiredPathParameter(name) {
    def value = params[name]?.toString()?.trim()
    if( !value ) {
        error "Missing required parameter --${name}"
    }
    return file(value, checkIfExists: true)
}

def barcodeSourceName(name, sampleId, modality) {
    if( !name.startsWith("${sampleId}.${modality}_") || !name.endsWith('.tsv') ) {
        return false
    }
    return !name.endsWith('_read_retention.tsv') &&
        !name.contains('.dual_tag_artifact_filter.')
}

def reportBarcodeSourceName(name) {
    return name.endsWith('.stats.tsv') ||
        name.endsWith('_sample_barcode.counts.tsv') ||
        name.endsWith('_barcode_gates.tsv') ||
        name.endsWith('_barcode_composition.tsv')
}

workflow {
    def scenario = params.scenario?.toString()?.trim()
    if( !(scenario in ['rna_only', 'dna_single', 'dna_dual']) ) {
        error "--scenario must be one of: rna_only, dna_single, dna_dual"
    }

    def repoRoot = requiredPathParameter('repo_root')
    def fixtureRoot = requiredPathParameter('fixture_root')
    def sourceRoot = requiredPathParameter('source_root')
    def resolvedOutdir = file(params.outdir).toAbsolutePath().normalize().toString()
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)

    def modality = scenario == 'rna_only' ? 'rna' : 'dna'
    def sampleId = scenario == 'rna_only' ? 'phase0_rna' : "phase0_${scenario}"
    def splitName = scenario == 'rna_only'
        ? 'phase0_rna_Normal'
        : "${sampleId}_Normal_H3K27ac"
    def tagmentation = scenario == 'rna_only' ? null : scenario.replace('dna_', '')
    def poisonPrefix = '/home/annan/micromamba/envs/tres/phase6-must-not-be-used'
    def taskMeta = [
        id                : "${modality}.${sampleId}.raw",
        tres_modality     : modality,
        tres_qc_stage     : 'raw_fastq',
        tres_split_name   : sampleId,
        runtime_env_prefix: poisonPrefix,
        runtime_tmpdir    : "${poisonPrefix}/tmp",
    ]

    def readsDir = new File(fixtureRoot.toString(), "reads/${scenario}")
    def readRoles = scenario == 'rna_only'
        ? ['I1.fastq', 'R1.fastq', 'R2.fastq']
        : scenario == 'dna_single'
            ? ['I1.fastq', 'I2.fastq', 'R1.fastq', 'R2.fastq']
            : ['I1.fastq', 'R1.fastq', 'R2.fastq']
    def rawReads = readRoles.collect { name ->
        file(new File(readsDir, name).canonicalPath, checkIfExists: true)
    }
    FASTQC(Channel.value(tuple(taskMeta, rawReads)))

    def bamMeta = taskMeta + [id: "${modality}.${splitName}.filtered_cells"]
    def bamInputs
    if( scenario == 'rna_only' ) {
        bamInputs = [tuple(
            bamMeta,
            file("${sourceRoot}/rna_align/${splitName}.filtered_cells.bam", checkIfExists: true),
            [],
            false
        )]
    }
    else {
        def alignedBam = requiredPathParameter('dna_aligned_bam')
        def alignedBai = requiredPathParameter('dna_aligned_bai')
        def dnaAlignDir = new File(sourceRoot.toString(), 'dna_align')
        bamInputs = [
            tuple(
                taskMeta + [id: "dna.${splitName}.aligned", tres_qc_stage: 'aligned', tres_split_name: splitName],
                alignedBam,
                alignedBai,
                true
            ),
            tuple(
                taskMeta + [id: "dna.${splitName}.markeddup", tres_qc_stage: 'markeddup', tres_split_name: splitName],
                file(new File(dnaAlignDir, "${splitName}_MarkedDup.bam").canonicalPath, checkIfExists: true),
                file(new File(dnaAlignDir, "${splitName}_MarkedDup.bam.bai").canonicalPath, checkIfExists: true),
                true
            ),
            tuple(
                taskMeta + [id: "dna.${splitName}.nodup", tres_qc_stage: 'nodup', tres_split_name: splitName],
                file(new File(dnaAlignDir, "${splitName}_NoDup.bam").canonicalPath, checkIfExists: true),
                file(new File(dnaAlignDir, "${splitName}_NoDup.bam.bai").canonicalPath, checkIfExists: true),
                true
            ),
        ]
    }
    def samtoolsQcScripts = [
        file("${repoRoot}/scripts/core_runtime/SamtoolsBamQc.sh", checkIfExists: true),
    ]
    SAMTOOLS_BAM_QC(Channel.fromList(bamInputs), samtoolsQcScripts)

    def statsDir = new File(sourceRoot.toString(), 'TrES_Stats')
    def barcodeFiles = statsDir.listFiles()
        .findAll { candidate -> candidate.isFile() && barcodeSourceName(candidate.name, sampleId, modality) }
        .sort { left, right -> left.name <=> right.name }
        .collect { candidate -> file(candidate.canonicalPath, checkIfExists: true) }

    def multiqcUpstream = new ArrayList(barcodeFiles)
    def reportUpstream = barcodeFiles
        .findAll { candidate -> reportBarcodeSourceName(candidate.name) }
        .collect()

    if( scenario == 'rna_only' ) {
        def starSummary = file("${sourceRoot}/rna_align/${splitName}.Solo.outGeneFull/Summary.csv", checkIfExists: true)
        def reportStarSummary = requiredPathParameter('rna_report_solo_summary')
        def starLog = file("${sourceRoot}/rna_align/${splitName}.Log.final.out", checkIfExists: true)
        multiqcUpstream.addAll([starSummary, starLog])
        reportUpstream.addAll([
            file("${sourceRoot}/TrES_Stats/${sampleId}.rna_read_retention.tsv", checkIfExists: true),
            file("${sourceRoot}/TrES_Stats/${splitName}.rna_filter_retention.tsv", checkIfExists: true),
            reportStarSummary,
            starLog,
            file("${sourceRoot}/pipeline_info/derived_contract/rna_sb_group_map.tsv", checkIfExists: true),
        ])
    }
    else {
        def duplicateMetrics = file("${sourceRoot}/dna_align/${splitName}.DuplicateMetrics.txt", checkIfExists: true)
        multiqcUpstream.add(duplicateMetrics)
        reportUpstream.addAll([
            file("${sourceRoot}/TrES_Stats/${sampleId}.dna_read_retention.tsv", checkIfExists: true),
            file("${sourceRoot}/TrES_Stats/${splitName}.dna_alignment_retention.tsv", checkIfExists: true),
            duplicateMetrics,
            file("${sourceRoot}/pipeline_info/derived_contract/dna_sb_group_map.tsv", checkIfExists: true),
            file("${sourceRoot}/pipeline_info/derived_contract/dna_mo_map.tsv", checkIfExists: true),
        ])
        if( scenario == 'dna_dual' ) {
            def artifactJson = file("${sourceRoot}/TrES_Stats/${sampleId}.dual_tag_artifact_filter.cutadapt.json", checkIfExists: true)
            def artifactSummary = file("${sourceRoot}/TrES_Stats/${sampleId}.dual_tag_artifact_filter.summary.tsv", checkIfExists: true)
            multiqcUpstream.addAll([artifactJson, artifactSummary])
            reportUpstream.add(artifactSummary)
        }
    }

    ch_multiqc_sources = Channel.fromList(multiqcUpstream)
        .mix(FASTQC.out.zip.map { meta, archive -> archive })
        .mix(SAMTOOLS_BAM_QC.out.flagstat.map { meta, report -> report })
        .mix(SAMTOOLS_BAM_QC.out.stats.map { meta, report -> report })
        .mix(SAMTOOLS_BAM_QC.out.idxstats.map { meta, report -> report })
        .mix(SAMTOOLS_BAM_QC.out.quickcheck.map { meta, report -> report })

    ch_multiqc_input = ch_multiqc_sources
        .collect()
        .map { files -> tuple(
            [id: 'tresflow'],
            files,
            file("${repoRoot}/assets/multiqc_config.yml", checkIfExists: true),
            [],
            [],
            []
        ) }

    def reportMeta = [
        id                       : 'tresflow',
        report_title             : scenario,
        pipeline_version         : 'v1.1.1',
        filter_dual_tag_artifacts: true,
        runtime_env_prefix       : poisonPrefix,
        runtime_tmpdir           : "${poisonPrefix}/tmp",
        samples                  : [[
            id              : sampleId,
            modality        : modality,
            dna_tagmentation: tagmentation,
            groups          : ['Normal'],
        ]],
    ]
    ch_report_input = Channel.fromList(reportUpstream)
        .mix(SAMTOOLS_BAM_QC.out.flagstat.map { meta, report -> report })
        .collect()
        .map { files -> tuple(reportMeta, files) }

    TRES_REPORT_HTML(
        ch_report_input,
        file("${repoRoot}/bin/render_tres_report.py", checkIfExists: true),
        file("${repoRoot}/lib/tresflow_qc", checkIfExists: true)
    )
    MULTIQC(ch_multiqc_input)
}
