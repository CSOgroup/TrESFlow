#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGN_DNA } from '../modules/local/align_dna/main'
include { FILTER_CANONICAL_DNA_ALIGNED_BAM } from '../modules/local/filter_canonical_dna_aligned_bam/main'
include { GATK4_MARKDUPLICATES } from '../modules/nf-core/gatk4/markduplicates/main'
include { NORMALIZE_DNA_MARKDUPLICATES } from '../modules/local/normalize_dna_markduplicates/main'
include { SPLIT_DUPLICATES_DNA } from '../modules/local/split_duplicates_dna/main'
include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main'

def nfcoreDnaMeta(splitName, meta, stage) {
    meta + [
        id             : splitName,
        tres_sample_id : meta.id,
        tres_split_name: splitName,
        tres_stage     : stage,
    ]
}

def restoreDnaMeta(meta) {
    meta + [id: meta.tres_sample_id ?: meta.id]
}

def hasMappedDnaReads(mappedReadsFile) {
    def value = mappedReadsFile.text.trim()
    if( !(value ==~ /[0-9]+/) ) {
        error "Invalid NoDup mapped-read count '${value}'"
    }
    value.toLong() > 0
}

workflow {
    def required = [
        'repo_root', 'outdir', 'bwa_reference', 'blacklist',
        'canonical_chromosomes', 'single_r1', 'single_r2', 'single_rg',
        'dual_r1', 'dual_r2', 'dual_rg', 'empty_markeddup_bam'
    ]
    for( name in required ) {
        if( !params[name] ) {
            error "Missing required parameter --${name}"
        }
    }

    def repoRoot = file(params.repo_root, checkIfExists: true)
    def resolvedOutdir = file(params.outdir).toAbsolutePath().normalize().toString()
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)

    def alignScripts = [
        file("${repoRoot}/scripts/core_runtime/AlignDNA.sh", checkIfExists: true),
    ]
    def canonicalScripts = [
        file("${repoRoot}/scripts/core_runtime/FilterCanonicalBam.sh", checkIfExists: true),
    ]
    def splitScripts = [
        file("${repoRoot}/scripts/core_runtime/SplitDuplicatesDNA.sh", checkIfExists: true),
    ]

    def bwaReference = params.bwa_reference as String
    def bwaReferenceName = new File(bwaReference).name
    def bwaIndexFiles = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac'].collect { suffix ->
        file("${bwaReference}${suffix}", checkIfExists: true)
    }
    def blacklist = file(params.blacklist, checkIfExists: true)
    def canonical = file(params.canonical_chromosomes, checkIfExists: true)
    def effectiveGenomeSize = '62000'
    def poisonedPrefix = '/home/annan/micromamba/envs/tres/phase5-must-not-be-used'

    def scenarios = [
        [
            sampleId: 'phase0_dna_single',
            splitName: 'phase0_dna_single_Normal_H3K27ac',
            r1: file(params.single_r1, checkIfExists: true),
            r2: file(params.single_r2, checkIfExists: true),
            rg: file(params.single_rg, checkIfExists: true),
        ],
        [
            sampleId: 'phase0_dna_dual',
            splitName: 'phase0_dna_dual_Normal_H3K27ac',
            r1: file(params.dual_r1, checkIfExists: true),
            r2: file(params.dual_r2, checkIfExists: true),
            rg: file(params.dual_rg, checkIfExists: true),
        ],
    ]

    ch_align = Channel.fromList(scenarios).map { scenario ->
        def meta = [
            id                       : scenario.sampleId,
            dna_effective_genome_size: effectiveGenomeSize,
            runtime_env_prefix       : poisonedPrefix,
            runtime_tmpdir           : "${poisonedPrefix}/tmp",
        ]
        tuple(
            scenario.splitName,
            meta,
            "${scenario.sampleId}_Normal",
            'H3K27ac',
            scenario.r1,
            scenario.r2,
            scenario.rg,
            bwaReferenceName,
            bwaIndexFiles,
            blacklist,
            effectiveGenomeSize
        )
    }

    ALIGN_DNA(ch_align, alignScripts)

    ch_canonical = ALIGN_DNA.out.bam
        .join(ALIGN_DNA.out.bai)
        .map { splitName, metaFromBam, bam, metaFromBai, bai ->
            tuple(splitName, metaFromBam, bam, bai, canonical)
        }
    FILTER_CANONICAL_DNA_ALIGNED_BAM(ch_canonical, canonicalScripts)

    ch_markdup = ALIGN_DNA.out.bam.map { splitName, meta, bam ->
        tuple(nfcoreDnaMeta(splitName, meta, 'markeddup'), bam)
    }
    GATK4_MARKDUPLICATES(ch_markdup, [], [])

    ch_normalize = GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai)
        .join(GATK4_MARKDUPLICATES.out.metrics)
        .map { nfMeta, bam, bai, metrics ->
            tuple(
                nfMeta.tres_split_name as String,
                restoreDnaMeta(nfMeta),
                bam,
                bai,
                metrics,
                canonical
            )
        }
    NORMALIZE_DNA_MARKDUPLICATES(ch_normalize, canonicalScripts)

    def emptyMeta = [
        id                       : 'phase5_empty',
        dna_effective_genome_size: effectiveGenomeSize,
        runtime_env_prefix       : poisonedPrefix,
        runtime_tmpdir           : "${poisonedPrefix}/tmp",
    ]
    ch_split = NORMALIZE_DNA_MARKDUPLICATES.out.bam
        .map { splitName, meta, bam ->
            tuple(splitName, meta, bam, effectiveGenomeSize)
        }
        .mix(Channel.value(tuple(
            'phase5_empty_Normal_H3K27ac',
            emptyMeta,
            file(params.empty_markeddup_bam, checkIfExists: true),
            effectiveGenomeSize
        )))
    SPLIT_DUPLICATES_DNA(ch_split, splitScripts)

    ch_coverage = SPLIT_DUPLICATES_DNA.out.bam
        .join(SPLIT_DUPLICATES_DNA.out.bai)
        .join(SPLIT_DUPLICATES_DNA.out.mapped_reads)
        .filter { splitName, metaFromBam, noDupBam, metaFromBai, noDupBai, metaFromCount, size, mappedReads ->
            hasMappedDnaReads(mappedReads)
        }
        .map { splitName, metaFromBam, noDupBam, metaFromBai, noDupBai, metaFromCount, size, mappedReads ->
            def coverageMeta = metaFromBam + [dna_effective_genome_size: size]
            tuple(nfcoreDnaMeta(splitName, coverageMeta, 'nodup_coverage'), noDupBam, noDupBai)
        }

    DEEPTOOLS_BAMCOVERAGE(ch_coverage, [], [], tuple([id: 'no_blacklist'], []))
}
