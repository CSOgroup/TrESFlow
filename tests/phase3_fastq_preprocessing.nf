#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIM_RNA_FASTQS } from '../modules/local/trim_rna_fastqs/main'
include { TRIM_DNA_FASTQS } from '../modules/local/trim_dna_fastqs/main'
include { DUAL_TAG_ARTIFACT_FILTER } from '../modules/local/dual_tag_artifact_filter/main'
include { COMPRESS_SPLIT_FASTQS as COMPRESS_RNA_SPLIT_FASTQS } from '../modules/local/compress_split_fastqs/main'
include { COMPRESS_SPLIT_FASTQS as COMPRESS_DNA_SPLIT_FASTQS } from '../modules/local/compress_split_fastqs/main'

workflow {
    if( !params.repo_root || !params.outdir ) {
        error 'Required parameters: --repo_root and --outdir'
    }

    def repoRoot = file(params.repo_root, checkIfExists: true)
    def fixtureRoot = file("${repoRoot}/tests/fixtures/fastq_preprocessing", checkIfExists: true)
    def resolvedOutdir = file(params.outdir).toAbsolutePath().normalize().toString()
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)

    def fixtureR1 = file("${fixtureRoot}/phase3_Normal_R1.fastq", checkIfExists: true)
    def fixtureR2 = file("${fixtureRoot}/phase3_Normal_R2.fastq", checkIfExists: true)
    def signatureFasta = file(
        "${repoRoot}/assets/dual_tag_artifact_23mers.fasta",
        checkIfExists: true
    )
    def trimHelper = file("${repoRoot}/bin/run_trim_galore.py", checkIfExists: true)
    def filterHelper = file(
        "${repoRoot}/bin/run_dual_tag_artifact_filter.py",
        checkIfExists: true
    )

    def poisonedPrefix = '/home/annan/micromamba/envs/tres/phase3-must-not-be-used'
    def commonMeta = [
        library_name      : 'PHASE3',
        runtime_env_prefix: poisonedPrefix,
        runtime_tmpdir    : "${poisonedPrefix}/tmp",
    ]

    def rnaMeta = commonMeta + [id: 'phase3_rna']
    TRIM_RNA_FASTQS(
        Channel.value(tuple('phase3_rna', rnaMeta, fixtureR1, fixtureR2)),
        trimHelper
    )

    def singleMeta = commonMeta + [id: 'phase3_dna_single', dna_tagmentation: 'single']
    def dualMeta = commonMeta + [id: 'phase3_dna_dual', dna_tagmentation: 'dual']
    TRIM_DNA_FASTQS(
        Channel.of(
            tuple('phase3_dna_single', singleMeta, fixtureR1, fixtureR2),
            tuple('phase3_dna_dual', dualMeta, fixtureR1, fixtureR2)
        ),
        trimHelper
    )

    ch_dual_trimmed = TRIM_DNA_FASTQS.out.trimmed.filter {
        _sampleId, meta, _trimmedR1, _trimmedR2 -> meta.dna_tagmentation == 'dual'
    }
    DUAL_TAG_ARTIFACT_FILTER(ch_dual_trimmed, signatureFasta, filterHelper)

    COMPRESS_RNA_SPLIT_FASTQS(
        Channel.value(tuple('phase3', rnaMeta, 'rna', [fixtureR1], [fixtureR2]))
    )
    COMPRESS_DNA_SPLIT_FASTQS(
        Channel.value(tuple('phase3', singleMeta, 'dna', [fixtureR1], [fixtureR2]))
    )
}
