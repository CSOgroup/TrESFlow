#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RNA_STARSOLO_ALIGN } from '../modules/local/rna_starsolo_align/main'
include { RNA_FILTERED_BAM } from '../modules/local/rna_filtered_bam/main'
include { RNA_COVERAGE } from '../modules/local/rna_coverage/main'

workflow {
    for( required in ['repo_root', 'outdir', 'usam', 'star_index', 'canonical_chromosomes', 'chrom_sizes'] ) {
        if( !params[required] ) {
            error "Missing required parameter --${required}"
        }
    }

    def repoRoot = file(params.repo_root, checkIfExists: true)
    def resolvedOutdir = file(params.outdir).toAbsolutePath().normalize().toString()
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)

    def starsoloScripts = [
        file("${repoRoot}/scripts/core_runtime/RNA_STARSOLO_ALIGN.sh", checkIfExists: true),
    ]
    def filteredBamScripts = [
        file("${repoRoot}/scripts/core_runtime/RNA_FILTERED_BAM.sh", checkIfExists: true),
        file("${repoRoot}/scripts/core_runtime/SummarizeRnaRetention.py", checkIfExists: true),
        file("${repoRoot}/scripts/core_runtime/FilterCanonicalBam.sh", checkIfExists: true),
    ]
    def coverageScripts = [
        file("${repoRoot}/scripts/core_runtime/RNA_COVERAGE.sh", checkIfExists: true),
    ]

    def splitName = 'phase0_rna_Normal'
    def poisonedPrefix = '/home/annan/micromamba/envs/tres/phase4-must-not-be-used'
    def meta = [
        id                : splitName,
        runtime_env_prefix: poisonedPrefix,
        runtime_tmpdir    : "${poisonedPrefix}/tmp",
    ]

    RNA_STARSOLO_ALIGN(
        Channel.value(tuple(
            splitName,
            meta,
            file(params.usam, checkIfExists: true),
            file(params.star_index, checkIfExists: true)
        )),
        starsoloScripts
    )

    ch_filtered_bam = RNA_STARSOLO_ALIGN.out.solo_dir
        .join(RNA_STARSOLO_ALIGN.out.aligned_bam)
        .map { name, soloMeta, soloDir, _bamMeta, alignedBam ->
            tuple(
                name,
                soloMeta,
                soloDir,
                alignedBam,
                file(params.canonical_chromosomes, checkIfExists: true)
            )
        }

    RNA_FILTERED_BAM(ch_filtered_bam, filteredBamScripts)

    ch_coverage = RNA_FILTERED_BAM.out.filtered_bam.map { name, filteredMeta, bam ->
        tuple(
            name,
            filteredMeta,
            bam,
            file(params.star_index, checkIfExists: true),
            file(params.chrom_sizes, checkIfExists: true)
        )
    }

    RNA_COVERAGE(ch_coverage, coverageScripts)
}
