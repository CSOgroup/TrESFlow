#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRESEQ } from './workflows/treseq'

def loadGroovySupportClass(sourcePath) {
    def loader = new groovy.lang.GroovyClassLoader()
    loader.parseClass(new java.io.File(sourcePath.toString()))
}

def parseNonNegativeIntegerParameter(rawValue, String parameterName) {
    def text = rawValue?.toString()?.trim()
    if( !text || !(text ==~ /[0-9]+/) ) {
        throw new IllegalArgumentException(
            "Parameter --${parameterName} must be a non-negative integer; received '${rawValue}'."
        )
    }
    def parsed = new BigInteger(text)
    if( parsed > Integer.MAX_VALUE ) {
        throw new IllegalArgumentException(
            "Parameter --${parameterName} exceeds the supported integer range: '${rawValue}'."
        )
    }
    return parsed.intValue()
}

def parseBooleanParameter(rawValue, String parameterName) {
    if( rawValue instanceof Boolean ) {
        return rawValue
    }
    def text = rawValue?.toString()?.trim()?.toLowerCase()
    if( text == 'true' ) {
        return true
    }
    if( text == 'false' ) {
        return false
    }
    throw new IllegalArgumentException(
        "Parameter --${parameterName} must be true or false; received '${rawValue}'."
    )
}

workflow {
    def runtimeSupport = loadGroovySupportClass("${projectDir}/lib/RuntimeSupport.groovy")
    def samplesheetParser = loadGroovySupportClass("${projectDir}/lib/SamplesheetParser.groovy")
    def workflowSupport = loadGroovySupportClass("${projectDir}/lib/WorkflowSupport.groovy")
    def rawSamplesheet = params.get('samplesheet')
    def rawOutdir = params.get('outdir')
    def rawCoreScriptsDir = params.get('core_scripts_dir')
    def avitiOpticalDuplicateDistance = null
    def filterDualTagArtifacts = null
    try {
        avitiOpticalDuplicateDistance = parseNonNegativeIntegerParameter(
            params.get('aviti_optical_duplicate_distance'),
            'aviti_optical_duplicate_distance'
        )
        filterDualTagArtifacts = parseBooleanParameter(
            params.get('filter_dual_tag_artifacts'),
            'filter_dual_tag_artifacts'
        )
    }
    catch( IllegalArgumentException e ) {
        error e.message
    }

    def resolvedSamplesheet = runtimeSupport.resolveLaunchPath(
        launchDir.toString(),
        rawSamplesheet
    )
    def resolvedOutdir = runtimeSupport.resolveLaunchPath(
        launchDir.toString(),
        rawOutdir ?: 'results'
    )
    def resolvedCoreScriptsDir = rawCoreScriptsDir
        ? runtimeSupport.resolveLaunchPath(launchDir.toString(), rawCoreScriptsDir)
        : runtimeSupport.resolveProjectPath(projectDir.toString(), 'scripts/core_runtime')
    def reportTitle = new File(resolvedOutdir).name
    def pipelineReleaseVersion = runtimeSupport.resolvePipelineReleaseVersion(
        projectDir.toString(),
        workflow.manifest.version
    )

    // Resolve launch-time paths once. Downstream modules consume these canonical
    // values, while repository-owned wrappers and assets continue to use projectDir.
    // Included modules have isolated params bindings, so run-scoped properties
    // expose the canonical output and core-script paths without changing channels.
    params.put('samplesheet', resolvedSamplesheet)
    params.put('outdir', resolvedOutdir)
    params.put('core_scripts_dir', resolvedCoreScriptsDir)
    params.put('aviti_optical_duplicate_distance', avitiOpticalDuplicateDistance)
    params.put('filter_dual_tag_artifacts', filterDualTagArtifacts)
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)
    java.lang.System.setProperty('tresflow.resolvedCoreScriptsDir', resolvedCoreScriptsDir)

    def deprecatedCliParams = [
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

    if( !resolvedSamplesheet ) {
        error "Missing required parameter: --samplesheet"
    }

    def samplesheetContract = null
    try {
        samplesheetContract = samplesheetParser.parseContract(
            resolvedSamplesheet,
            [
                outdir          : resolvedOutdir,
                barcode_defaults: params.barcode_defaults,
            ]
        )
    }
    catch( IllegalArgumentException e ) {
        error e.message
    }

    def runtimeConfig = samplesheetContract['runtime'] as Map
    def referenceConfig = samplesheetContract['references'] as Map
    def modalityConfig = samplesheetContract['modalities'] as Map
    def runtimeParams = [
        runtime_env_prefix: runtimeConfig['env_prefix'],
        runtime_tmpdir    : runtimeConfig['tmpdir'],
    ]
    def sampleRows = samplesheetContract['samples'] as List<Map>

    log.warn """
    ================================================================================
    TrESFlow runtime TMPDIR resolved for this run:
      ${runtimeParams.runtime_tmpdir}

    This directory can become very large on production FASTQ/BAM runs. Monitor free
    space on the filesystem that backs this path.
    ================================================================================
    """.stripIndent().trim()

    runtimeSupport.validateRuntimeContract(runtimeParams)
    runtimeSupport.validateConfiguredDirectory('core scripts dir', resolvedCoreScriptsDir)
    def codonPreflightOutput = runtimeSupport.runCodonSeqPreflight(
        runtimeParams,
        projectDir.toString()
    )
    workflowSupport.validateReferenceContract(
        referenceConfig,
        modalityConfig,
        sampleRows
    )
    def canonicalChromosomeContracts = runtimeSupport.writeCanonicalChromosomeContracts(
        runtimeParams,
        projectDir.toString(),
        resolvedOutdir,
        referenceConfig,
        modalityConfig
    )
    sampleRows.each { row ->
        def chromosomeContract = canonicalChromosomeContracts[row.modality]
        if( !chromosomeContract ) {
            error "Missing canonical chromosome contract for modality '${row.modality}'"
        }
        row.canonical_chromosomes = chromosomeContract.allowlist
        row.canonical_chrom_sizes = chromosomeContract.chrom_sizes
        row.chromosome_naming = chromosomeContract.style
    }
    canonicalChromosomeContracts.each { modality, contract ->
        log.info "Resolved ${modality.toUpperCase()} canonical chromosomes " +
            "(${contract.style}): ${contract.contigs.join(', ')}"
    }
    runtimeSupport.writeRuntimeContract(
        resolvedOutdir,
        runtimeSupport.configuredRuntimeTools(runtimeParams),
        codonPreflightOutput,
        runtimeSupport.runtimeContext(runtimeParams)
    )

    TRESEQ(
        sampleRows,
        [
            report_title    : reportTitle,
            pipeline_version: pipelineReleaseVersion,
        ]
    )
}
