#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BARCODE_GATE_METRICS } from '../modules/local/barcode_gate_metrics/main'
include { TRES_REPORT_HTML } from '../modules/local/tres_report_html/main'

workflow {
    if( !params.repo_root || !params.fixture_dir || !params.outdir ) {
        error 'Required parameters: --repo_root, --fixture_dir, and --outdir'
    }

    def repoRoot = file(params.repo_root, checkIfExists: true)
    def fixtureRoot = file(params.fixture_dir, checkIfExists: true)
    def resolvedOutdir = file(params.outdir).toAbsolutePath().normalize().toString()
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)

    def r1 = file("${fixtureRoot}/R1.fastq", checkIfExists: true)
    def r2 = file("${fixtureRoot}/R2.fastq", checkIfExists: true)
    def tagRecords = file("${fixtureRoot}/tag_records.tsv", checkIfExists: true)
    def sbGroupMap = file("${fixtureRoot}/dna_sb_group_map.tsv", checkIfExists: true)
    def moMap = file("${fixtureRoot}/dna_mo_map.tsv", checkIfExists: true)
    def deliberatelyInvalidRuntime = '/home/annan/micromamba/envs/tres/phase1-must-not-be-used'
    def barcodeMeta = [
        id: 'sample',
        runtime_env_prefix: deliberatelyInvalidRuntime,
        runtime_tmpdir: '/home/annan/micromamba/envs/tres/phase1-must-not-be-used/tmp',
    ]
    def barcodeHelpers = [
        file("${repoRoot}/bin/write_barcode_gate_metrics.py", checkIfExists: true),
        file("${repoRoot}/bin/run_split_reads_dna.py", checkIfExists: true),
        file("${repoRoot}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]

    BARCODE_GATE_METRICS(
        Channel.value(tuple('sample', barcodeMeta, 'dna', r1, r2, tagRecords, [sbGroupMap, moMap])),
        barcodeHelpers
    )

    def reportStaticInputs = file("${fixtureRoot}/report_static/*", checkIfExists: true)
    def reportMeta = [
        id: 'phase1-python-helpers',
        report_title: 'Phase 1 portable helpers',
        pipeline_version: 'v1.1.1',
        filter_dual_tag_artifacts: false,
        runtime_env_prefix: deliberatelyInvalidRuntime,
        runtime_tmpdir: '/home/annan/micromamba/envs/tres/phase1-must-not-be-used/tmp',
        samples: [[
            id: 'sample',
            modality: 'dna',
            dna_tagmentation: 'single',
            groups: ['run-alpha', 'run-beta'],
        ]],
    ]
    ch_report_input = BARCODE_GATE_METRICS.out.metrics.map {
        sampleId, _meta, gates, composition ->
        tuple(reportMeta, [gates, composition, sbGroupMap, moMap] + reportStaticInputs)
    }

    TRES_REPORT_HTML(
        ch_report_input,
        file("${repoRoot}/bin/render_tres_report.py", checkIfExists: true),
        file("${repoRoot}/lib/tresflow_qc", checkIfExists: true)
    )
}
