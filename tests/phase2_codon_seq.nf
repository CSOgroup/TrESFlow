#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TAG_RNA_SAMPLE_BARCODE } from '../modules/local/tag_rna_sb/main'
include { TAG_RNA_UMI } from '../modules/local/tag_rna_umi/main'
include { TAG_RNA_CELL_BARCODE } from '../modules/local/tag_rna_cell_barcode/main'
include { SPLIT_RNA_READS } from '../modules/local/split_rna_reads/main'
include { FQ_TO_SAM } from '../modules/local/fq_to_sam/main'
include { TAG_DNA_SAMPLE_BARCODE } from '../modules/local/tag_dna_sb/main'
include { TAG_DNA_MODALITY_BARCODE } from '../modules/local/tag_dna_modality/main'
include { TAG_DNA_CELL_BARCODE } from '../modules/local/tag_dna_cell_barcode/main'
include { SPLIT_DNA_READS } from '../modules/local/split_dna_reads/main'

def asPathList(value) {
    value instanceof List ? value : [value]
}

workflow {
    if( !params.repo_root || !params.outdir ) {
        error 'Required parameters: --repo_root and --outdir'
    }

    def repoRoot = file(params.repo_root, checkIfExists: true)
    def fixtureRoot = file("${repoRoot}/tests/fixtures/codon_seq", checkIfExists: true)
    def testdataRoot = file("${repoRoot}/assets/testdata", checkIfExists: true)
    def coreRoot = file("${repoRoot}/scripts/core_runtime", checkIfExists: true)
    def resolvedOutdir = file(params.outdir).toAbsolutePath().normalize().toString()
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)

    def invalidHostRuntime = '/home/annan/micromamba/envs/tres/phase2-must-not-be-used'
    def commonMeta = [
        library_name: 'PHASE2',
        runtime_env_prefix: invalidHostRuntime,
        runtime_tmpdir: "${invalidHostRuntime}/tmp",
    ]

    def tagHelpers = [
        file("${repoRoot}/bin/run_tag.py", checkIfExists: true),
        file("${repoRoot}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def tagCodon = [
        file("${coreRoot}/Tag.codon", checkIfExists: true),
        file("${coreRoot}/utils.codon", checkIfExists: true),
    ]
    def umiHelpers = [
        file("${repoRoot}/bin/run_tag_umi.py", checkIfExists: true),
        file("${repoRoot}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def umiCodon = [file("${coreRoot}/Tag_UMI.codon", checkIfExists: true)]
    def cellHelpers = [
        file("${repoRoot}/bin/run_tag_lig3.py", checkIfExists: true),
        file("${repoRoot}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def cellCodon = [
        file("${coreRoot}/Tag_Lig3.codon", checkIfExists: true),
        file("${coreRoot}/utils.codon", checkIfExists: true),
    ]
    def rnaSplitHelpers = [
        file("${repoRoot}/bin/run_split_reads_rna.py", checkIfExists: true),
        file("${repoRoot}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def dnaSplitHelpers = [
        file("${repoRoot}/bin/run_split_reads_dna.py", checkIfExists: true),
        file("${repoRoot}/bin/tresflow_fastq_utils.py", checkIfExists: true),
    ]
    def splitCodon = [file("${coreRoot}/Split_ReadsV2.codon", checkIfExists: true)]
    def fqToSamHelpers = [file("${repoRoot}/bin/run_fq_to_sam.py", checkIfExists: true)]
    def fqToSamCodon = [file("${coreRoot}/FqToSAM.codon", checkIfExists: true)]

    def cellWhitelist = file("${testdataRoot}/TrESFlow_References/ligation_barcode_whitelist.txt", checkIfExists: true)
    def rnaSbMap = file("${fixtureRoot}/rna_sb_group_map.tsv", checkIfExists: true)
    def rnaI1 = file("${testdataRoot}/test_rna_I1.fastq", checkIfExists: true)
    def rnaR1 = file("${testdataRoot}/test_rna_R1.fastq", checkIfExists: true)
    def rnaR2 = file("${testdataRoot}/test_rna_R2.fastq", checkIfExists: true)
    def rnaMeta = commonMeta + [
        id: 'phase2_rna',
        rna_sb_barcode_source: 'r2',
        rna_sb_barcode_len: 4,
        sample_bc_len: 4,
        sample_bc_start: 0,
        sample_hd: 0,
        sample_tag: 'SB',
        sample_first_pass: 'first_pass',
        sample_reverse_complement: 'rev',
        umi_bc_len: 10,
        umi_bc_start: 4,
        umi_tag: 'UM',
        cell_bc_len: 8,
        cell_hd: 1,
        cell_tag: 'CB',
    ]

    TAG_RNA_SAMPLE_BARCODE(
        Channel.value(tuple('phase2_rna', rnaMeta, rnaR1, rnaR2, rnaSbMap)),
        tagHelpers,
        tagCodon
    )
    ch_rna_umi = TAG_RNA_SAMPLE_BARCODE.out.tagged.map {
        sampleId, meta, taggedR1, taggedR2, readSetCounts ->
        tuple(sampleId, meta, rnaR2, taggedR1, taggedR2, readSetCounts)
    }
    TAG_RNA_UMI(ch_rna_umi, umiHelpers, umiCodon)
    ch_rna_cell = TAG_RNA_UMI.out.tagged.map {
        sampleId, meta, taggedR1, taggedR2, readSetCounts ->
        tuple(sampleId, meta, rnaI1, taggedR1, taggedR2, cellWhitelist, readSetCounts)
    }
    TAG_RNA_CELL_BARCODE(ch_rna_cell, cellHelpers, cellCodon)
    ch_rna_split = TAG_RNA_CELL_BARCODE.out.tagged.map { sampleId, meta, taggedR1, taggedR2 ->
        tuple(sampleId, meta, taggedR1, taggedR2, rnaSbMap)
    }
    SPLIT_RNA_READS(ch_rna_split, rnaSplitHelpers, splitCodon)
    ch_fq_to_sam = SPLIT_RNA_READS.out.split_fastqs.flatMap { sampleId, meta, splitR1s, splitR2s ->
        def r2ByName = asPathList(splitR2s).collectEntries { path ->
            [(path.getName().replaceFirst('_R2\\.fastq$', '')): path]
        }
        asPathList(splitR1s).collect { splitR1 ->
            def splitName = splitR1.getName().replaceFirst('_R1\\.fastq$', '')
            tuple(splitName, meta, splitR1, r2ByName[splitName])
        }
    }
    FQ_TO_SAM(ch_fq_to_sam, fqToSamHelpers, fqToSamCodon)

    def dnaSbMap = file("${fixtureRoot}/dna_sb_group_map.tsv", checkIfExists: true)
    def dnaMoMap = file("${fixtureRoot}/dna_mo_map.tsv", checkIfExists: true)
    def dnaR1 = file("${testdataRoot}/test_dna_R1.fastq", checkIfExists: true)
    def dnaR2 = file("${testdataRoot}/test_dna_R2.fastq", checkIfExists: true)
    def dnaSingleI1 = file("${testdataRoot}/test_dna_I1.fastq", checkIfExists: true)
    def dnaSingleI2 = file("${testdataRoot}/test_dna_I2_single_lig.fastq", checkIfExists: true)
    def dnaDualI1 = file("${testdataRoot}/test_dna_I1_dual_lig.fastq", checkIfExists: true)
    def dnaSingleWhitelist = file("${fixtureRoot}/dna_single_modality_whitelist.txt", checkIfExists: true)
    def dnaDualWhitelist = file("${fixtureRoot}/dna_dual_modality_whitelist.txt", checkIfExists: true)
    def dnaCommon = commonMeta + [
        sample_bc_len: 4,
        sample_hd: 0,
        sample_tag: 'SB',
        sample_first_pass: 'first_pass',
        sample_reverse_complement: 'rev',
        modality_bc_len: 8,
        modality_hd: 1,
        modality_tag: 'MO',
        modality_first_pass: 'not_first_pass',
        modality_reverse_complement: 'rev',
        cell_bc_len: 8,
        cell_hd: 1,
        cell_tag: 'CB',
        dna_ligation_index_read: 'i1',
    ]
    def dnaSingleMeta = dnaCommon + [
        id: 'phase2_dna_single',
        dna_tagmentation: 'single',
        dna_sb_barcode_source: 'i2',
        dna_sb_barcode_len: 4,
        dna_sample_index_read: 'i2',
        dna_modality_index_read: 'i2',
        dna_ligation_start_positions: '15,53,91',
        sample_bc_start: 14,
        modality_bc_start: 18,
    ]
    def dnaDualMeta = dnaCommon + [
        id: 'phase2_dna_dual',
        dna_tagmentation: 'dual',
        dna_sb_barcode_source: 'i1',
        dna_sb_barcode_len: 3,
        dna_sample_index_read: 'i1',
        dna_modality_index_read: 'i1',
        dna_ligation_start_positions: '41,79,117',
        sample_bc_len: 3,
        sample_bc_start: 0,
        sample_reverse_complement: 'fw',
        modality_bc_start: 3,
        modality_reverse_complement: 'fw',
    ]

    ch_dna_sb = Channel.of(
        tuple('phase2_dna_single', dnaSingleMeta, dnaSingleI2, dnaR1, dnaR2, dnaSbMap),
        tuple('phase2_dna_dual', dnaDualMeta, dnaDualI1, dnaR1, dnaR2, dnaSbMap)
    )
    TAG_DNA_SAMPLE_BARCODE(ch_dna_sb, tagHelpers, tagCodon)

    ch_dna_mo_context = Channel.of(
        tuple('phase2_dna_single', dnaSingleMeta, dnaSingleI2, dnaSingleWhitelist),
        tuple('phase2_dna_dual', dnaDualMeta, dnaDualI1, dnaDualWhitelist)
    )
    ch_dna_mo = ch_dna_mo_context.join(TAG_DNA_SAMPLE_BARCODE.out.tagged).map {
        sampleId, meta, indexRead, modalityWhitelist, taggedMeta, taggedR1, taggedR2, readSetCounts ->
        tuple(sampleId, meta, indexRead, taggedR1, taggedR2, modalityWhitelist, readSetCounts)
    }
    TAG_DNA_MODALITY_BARCODE(ch_dna_mo, tagHelpers, tagCodon)

    ch_dna_cell_context = Channel.of(
        tuple('phase2_dna_single', dnaSingleMeta, dnaSingleI1),
        tuple('phase2_dna_dual', dnaDualMeta, dnaDualI1)
    )
    ch_dna_cell = ch_dna_cell_context.join(TAG_DNA_MODALITY_BARCODE.out.tagged).map {
        sampleId, meta, ligationRead, taggedMeta, taggedR1, taggedR2, readSetCounts ->
        tuple(sampleId, meta, ligationRead, taggedR1, taggedR2, cellWhitelist, readSetCounts)
    }
    TAG_DNA_CELL_BARCODE(ch_dna_cell, cellHelpers, cellCodon)
    ch_dna_split = TAG_DNA_CELL_BARCODE.out.tagged.map { sampleId, meta, taggedR1, taggedR2 ->
        tuple(sampleId, meta, taggedR1, taggedR2, dnaMoMap, dnaSbMap)
    }
    SPLIT_DNA_READS(ch_dna_split, dnaSplitHelpers, splitCodon)
}
