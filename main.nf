#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import RuntimeSupport
import SamplesheetParser
import WorkflowSupport

final String resolvedOutdir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.outdir ?: 'results')
final String resolvedCoreScriptsDir = RuntimeSupport.resolveProjectPath(projectDir.toString(), params.core_scripts_dir ?: 'scripts/core_runtime')

final Map<String, String> deprecatedCliParams = [
    runtime_env_prefix          : 'runtime.env_prefix',
    runtime_tmpdir              : 'runtime.tmpdir',
    ligation_barcode_whitelist : 'references.ligation_barcode_whitelist',
    rna_ref_base_dir           : 'references.rna_ref_dir',
    rna_align_species          : 'references.species',
    rna_ref_dir                : 'references.rna_ref_dir',
    dna_ref_dir                : 'references.dna_ref_dir',
    dna_bwa_reference          : 'the inferred prefix from references.dna_ref_dir',
    dna_blacklist_bed          : 'references.dna_blacklist_bed',
    dna_chrom_sizes            : 'references.dna_chrom_sizes',
    dna_effective_genome_size  : 'references.dna_effective_genome_size',
]

deprecatedCliParams.each { paramName, replacement ->
    if( params.containsKey(paramName) && params[paramName]?.toString()?.trim() ) {
        error "Deprecated parameter --${paramName} is no longer supported. Configure ${replacement} in the samplesheet instead."
    }
}

if( !params.samplesheet ) {
    error "Missing required parameter: --samplesheet"
}

Map samplesheetContract = null
try {
    samplesheetContract = SamplesheetParser.parseContract(
        params.samplesheet as String,
        [
            outdir          : resolvedOutdir,
            barcode_defaults: params.barcode_defaults,
        ]
    )
}
catch( IllegalArgumentException e ) {
    error e.message
}

final Map runtimeConfig = samplesheetContract['runtime'] as Map
final Map referenceConfig = samplesheetContract['references'] as Map
final Map modalityConfig = samplesheetContract['modalities'] as Map
final Map runtimeParams = [
    runtime_env_prefix: runtimeConfig['env_prefix'],
    runtime_tmpdir    : runtimeConfig['tmpdir'],
]
final List<Map> sampleRows = samplesheetContract['samples'] as List<Map>

// The runtime contract is enforced once up front so every downstream task sees
// the same validated toolchain and the same pinned Codon/Seq preflight result.
def runCodonSeqPreflight(final Map runtimeParams) {
    final File preflight = new File(projectDir.toString(), 'bin/check_codon_seq_host.sh')
    final String codonBin = RuntimeSupport.runtimeToolPath(runtimeParams, 'codon')
    final String codonHome = RuntimeSupport.runtimeCodonHome(runtimeParams)

    if( !preflight.exists() ) {
        throw new IllegalStateException(
            "Global pinned Codon/Seq preflight failed. Missing required preflight script: ${preflight}"
        )
    }

    final ProcessBuilder processBuilder = new ProcessBuilder('bash', preflight.toString())
        .directory(projectDir.toFile())
        .redirectErrorStream(true)

    final Map<String, String> env = processBuilder.environment()
    if( codonBin ) {
        env.put('CODON_BIN', codonBin)
    }
    if( codonHome ) {
        env.put('CODON_HOME', codonHome)
    }

    final Process process = processBuilder.start()

    final String output = process.inputStream.getText('UTF-8').trim()
    final int exitCode = process.waitFor()

    if( exitCode != 0 ) {
        final String detail = output ? "\n${output}" : ''
        throw new IllegalStateException(
            "Global pinned Codon/Seq preflight failed. Every pipeline run requires codon 0.16.3 and Seq 0.11.3.${detail}"
        )
    }

    return output
}

log.warn """
================================================================================
TrESFlow runtime TMPDIR resolved for this run:
  ${runtimeParams.runtime_tmpdir}

This directory can become very large on production FASTQ/BAM runs. Monitor free
space on the filesystem that backs this path.
================================================================================
""".stripIndent().trim()

RuntimeSupport.validateRuntimeContract(runtimeParams)
RuntimeSupport.validateConfiguredDirectory('core scripts dir', resolvedCoreScriptsDir)
final String codonPreflightOutput = runCodonSeqPreflight(runtimeParams)
WorkflowSupport.validateReferenceContract(
    referenceConfig,
    modalityConfig,
    sampleRows
)
RuntimeSupport.writeRuntimeContract(
    resolvedOutdir,
    RuntimeSupport.configuredRuntimeTools(runtimeParams),
    codonPreflightOutput,
    RuntimeSupport.runtimeContext(runtimeParams)
)

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ(sampleRows)
}
