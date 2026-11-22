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
    if( params.help ) {
        exit(0, '''
        Usage: nextflow run . --samplesheet <yaml> --outdir <directory> [options]

        Input modes:
          Raw reads from the samplesheet                 Run the full pipeline (default).
          --bulk true                                    Count barcoded bulk RNA and/or DNA targets.
          --bulk_reference_bed <path>                    Optional BED for bulk DNA fragment overlap QC.
          --preprocess_only true                         Publish reference-free preprocessing outputs.
          --preprocessed_input true                      Analyze split FASTQ pairs listed in samples.<sample>.preprocessed.
          --preprocessed_input <directory|manifest.json> Restart from a version 1 manifest.

        FASTQ paths are relative to the samplesheet directory. Supply the original
        sample and group labels, a mark for DNA, and the runtime and genome references.
        Set mark_barcode on a DNA split to check its MO tags.
        Use a separate output directory for analysis.
        Use either --preprocess_only true or --preprocessed_input.
        '''.stripIndent().trim())
    }
    if( params.containsKey('restart_from_preprocessed') ) {
        error 'Use --preprocessed_input in place of --restart_from_preprocessed.'
    }

    def runtimeSupport = loadGroovySupportClass("${projectDir}/lib/RuntimeSupport.groovy")
    def samplesheetParser = loadGroovySupportClass("${projectDir}/lib/SamplesheetParser.groovy")
    def workflowSupport = loadGroovySupportClass("${projectDir}/lib/WorkflowSupport.groovy")
    def rawSamplesheet = params.get('samplesheet')
    def rawPreprocessedInput = params.get('preprocessed_input')
    def rawOutdir = params.get('outdir')
    def rawCoreScriptsDir = params.get('core_scripts_dir')
    def avitiOpticalDuplicateDistance = null
    def filterDualTagArtifacts = null
    def bulk = null
    def preprocessOnly = null
    try {
        bulk = parseBooleanParameter(params.bulk, 'bulk')
        preprocessOnly = parseBooleanParameter(params.get('preprocess_only'), 'preprocess_only')
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
    if( rawPreprocessedInput != null && !(rawPreprocessedInput instanceof Boolean) &&
        (!(rawPreprocessedInput instanceof String) || !rawPreprocessedInput.trim()) ) {
        error 'Parameter --preprocessed_input accepts true, false, a preprocessing output directory, or a version 1 restart manifest path.'
    }
    def inputMode = rawPreprocessedInput?.toString()?.trim()?.toLowerCase()
    def fastqRestart = inputMode == 'true'
    def manifestPath = inputMode && !(inputMode in ['true', 'false'])
        ? runtimeSupport.resolveLaunchPath(launchDir.toString(), rawPreprocessedInput)
        : null
    if( (fastqRestart || manifestPath) && preprocessOnly ) {
        error 'Use either --preprocessed_input or --preprocess_only true.'
    }
    if( bulk && (preprocessOnly || manifestPath) ) {
        error 'Use --bulk true with raw reads or --preprocessed_input true.'
    }
    def referenceBed = params.bulk_reference_bed ? runtimeSupport.resolveLaunchPath(launchDir.toString(), params.bulk_reference_bed) : null
    if( referenceBed && (!bulk || (!new File(referenceBed).isFile() || !new File(referenceBed).canRead())) ) {
        error '--bulk_reference_bed requires --bulk true and a readable BED file.'
    }
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
    params.put('preprocess_only', preprocessOnly)
    if( preprocessOnly ) {
        log.info 'Preprocessing-only mode: enabling gzip split FASTQ publication; --publish_split_fastqs false is overridden. Alignment and downstream BAM/coverage stages will not run.'
    }
    java.lang.System.setProperty('tresflow.resolvedOutdir', resolvedOutdir)
    java.lang.System.setProperty('tresflow.resolvedWorkDir', workflow.workDir.toString())
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
    def restart = null
    try {
        // Check source paths before the samplesheet parser writes derived files.
        def manifestSupport = (fastqRestart || manifestPath)
            ? loadGroovySupportClass("${projectDir}/lib/PreprocessedManifest.groovy")
            : null
        if( manifestPath ) restart = manifestSupport.read(manifestPath, resolvedOutdir, workflow.workDir.toString())
        samplesheetContract = samplesheetParser.parseContract(
            resolvedSamplesheet,
            [
                preprocess_only : preprocessOnly,
                restart         : !!manifestPath,
                fastq_restart   : fastqRestart,
                outdir          : resolvedOutdir,
                barcode_defaults: params.barcode_defaults,
                sb_oligo_lookup : params.sb_oligo_lookup ?: "${projectDir}/assets/sb_oligo_lookup.v1.tsv",
            ]
        )
        if( fastqRestart ) {
            restart = manifestSupport.readFastqList(
                samplesheetContract.splits, resolvedSamplesheet, resolvedOutdir, workflow.workDir.toString()
            )
        }
        else if( restart ) manifestSupport.bindSamples(restart, samplesheetContract.samples)
    }
    catch( IllegalArgumentException e ) {
        error e.message
    }

    // Emit validated identity declarations on every launch, including resume,
    // before runtime/reference preflight or any process submission.
    samplesheetContract.sb_identity_warnings.each { warning ->
        log.warn samplesheetParser.delimitIdentityWarning(warning)
    }

    def runtimeConfig = samplesheetContract['runtime'] as Map
    def referenceConfig = samplesheetContract['references'] as Map
    def modalityConfig = samplesheetContract['modalities'] as Map
    def runtimeParams = [
        runtime_env_prefix: runtimeConfig['env_prefix'],
        runtime_tmpdir    : runtimeConfig['tmpdir'],
        preprocess_only   : preprocessOnly,
        bulk_dna          : bulk && modalityConfig.dna,
    ]
    def sampleRows = samplesheetContract['samples'] as List<Map>

    sampleRows.each { row ->
        row.bulk = bulk
        row.bulk_reference_bed = referenceBed
    }
    if( referenceBed && !modalityConfig.dna ) error '--bulk_reference_bed requires DNA inputs.'
    if( bulk && modalityConfig.rna ) {
        ['exonGeTrInfo.tab', 'geneInfo.tab', 'chrStart.txt'].each { filename ->
            if( !new File(referenceConfig.rna_ref_dir.toString(), filename).isFile() ) {
                error "Bulk RNA requires an annotated STAR index containing ${filename}: ${referenceConfig.rna_ref_dir}"
            }
        }
    }

    log.warn """
    TrESFlow storage paths:
      results: ${resolvedOutdir}
      TMPDIR:  ${runtimeParams.runtime_tmpdir}
      workDir: ${workflow.workDir}

    Large runs can use substantial disk space; monitor free space.
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
        sampleRows,
        preprocessOnly,
        !!restart
    )
    if( !preprocessOnly ) {
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
            preprocess_only : preprocessOnly,
            bulk            : bulk,
            report_title    : reportTitle,
            pipeline_version: pipelineReleaseVersion,
            restart         : restart,
        ]
    )
}
