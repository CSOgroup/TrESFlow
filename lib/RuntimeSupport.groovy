import groovy.json.JsonSlurper

class RuntimeSupport {

    private static final List<Map> STANDARD_RUNTIME_TOOLS = [
        [name: 'python3', binary: 'python3'],
        [name: 'samtools', binary: 'samtools'],
        [name: 'bwa-mem2', binary: 'bwa-mem2'],
        [name: 'bamCoverage', binary: 'bamCoverage'],
        [name: 'gatk', binary: 'gatk'],
    ]

    static void validateConfiguredExecutable(final String label, final String rawPath) {
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            throw new IllegalStateException("Missing configured executable path for ${label}")
        }

        final File executable = new File(path)
        if( !executable.exists() || !executable.canExecute() ) {
            throw new IllegalStateException(
                "Configured executable for ${label} is missing or not executable: ${executable}"
            )
        }
    }

    static void validateConfiguredDirectory(final String label, final String rawPath) {
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            throw new IllegalStateException("Missing configured directory path for ${label}")
        }

        final File directory = new File(path)
        if( !directory.exists() || !directory.isDirectory() ) {
            throw new IllegalStateException(
                "Configured directory for ${label} is missing or not a directory: ${directory}"
            )
        }
    }

    static String runtimeEnvPrefix(final Map params) {
        return (params.runtime_env_prefix ?: '').toString().trim()
    }

    static String runtimeTmpdir(final Map params) {
        return (params.runtime_tmpdir ?: '').toString().trim()
    }

    static String runtimeBinDir(final Map params) {
        final String envPrefix = runtimeEnvPrefix(params)
        return envPrefix ? "${envPrefix}/bin" : ''
    }

    static String runtimeToolPath(final Map params, final String binary) {
        final String binDir = runtimeBinDir(params)
        return binDir ? "${binDir}/${binary}" : ''
    }

    static String resolvePath(final String rawBaseDir, final Object rawPath) {
        final String baseDir = rawBaseDir?.toString()?.trim()
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            return path
        }

        final File candidate = new File(path)
        if( candidate.isAbsolute() || !baseDir ) {
            return candidate.canonicalPath
        }

        return new File(baseDir, path).canonicalPath
    }

    static String resolveLaunchPath(final String rawLaunchDir, final Object rawPath) {
        return resolvePath(rawLaunchDir, rawPath)
    }

    static String resolveProjectPath(final String rawProjectDir, final Object rawPath) {
        return resolvePath(rawProjectDir, rawPath)
    }

    static String resolvePipelineReleaseVersion(
        final String rawProjectDir,
        final Object manifestVersion
    ) {
        final File projectDirectory = new File(rawProjectDir).canonicalFile
        final File resolver = new File(projectDirectory, 'bin/resolve_tresflow_release_version.sh')
        if( !resolver.exists() ) {
            throw new IllegalStateException(
                "Missing repository release-version resolver: ${resolver}"
            )
        }

        final Process process = new ProcessBuilder(
            'bash',
            resolver.canonicalPath,
            projectDirectory.canonicalPath,
            (manifestVersion ?: '').toString()
        )
            .directory(projectDirectory)
            .redirectErrorStream(true)
            .start()
        final String output = process.inputStream.getText('UTF-8').trim()
        final int exitCode = process.waitFor()
        if( exitCode != 0 || !output ) {
            throw new IllegalStateException(
                "Unable to resolve the TrESFlow release version: ${output ?: 'no output'}"
            )
        }
        return output
    }

    static Map writeCanonicalChromosomeContracts(
        final Map runtimeParams,
        final String rawProjectDir,
        final String rawOutdir,
        final Map references,
        final Map modalities
    ) {
        final File projectDirectory = new File(rawProjectDir).canonicalFile
        final File resolver = new File(projectDirectory, 'bin/resolve_canonical_chromosomes.py')
        final File outputDirectory = new File(
            new File(rawOutdir).canonicalFile,
            'pipeline_info/derived_contract'
        )
        final String pythonBin = runtimeToolPath(runtimeParams, 'python3')

        validateConfiguredExecutable('canonical chromosome resolver', resolver.canonicalPath)
        validateConfiguredExecutable('runtime python3', pythonBin)

        final List<String> command = [
            pythonBin,
            resolver.canonicalPath,
            '--output-dir',
            outputDirectory.canonicalPath,
        ]

        if( modalities.rna as boolean ) {
            command.addAll([
                '--rna-chrom-sizes',
                references.rna_chrom_sizes.toString(),
            ])
        }
        if( modalities.dna as boolean ) {
            command.addAll([
                '--dna-bwa-ann',
                "${references.dna_bwa_reference}.ann".toString(),
            ])
            final String dnaChromSizes = references.dna_chrom_sizes?.toString()?.trim()
            if( dnaChromSizes ) {
                command.addAll(['--dna-chrom-sizes', dnaChromSizes])
            }
        }

        final Process process = new ProcessBuilder(command)
            .directory(projectDirectory)
            .redirectErrorStream(true)
            .start()
        final String output = process.inputStream.getText('UTF-8').trim()
        final int exitCode = process.waitFor()
        if( exitCode != 0 ) {
            throw new IllegalArgumentException(
                "Canonical chromosome resolution failed for the configured reference index: ${output}"
            )
        }

        try {
            return new JsonSlurper().parseText(output) as Map
        }
        catch( Exception error ) {
            throw new IllegalStateException(
                "Canonical chromosome resolver returned invalid output: ${output}",
                error
            )
        }
    }

    static List<Map> standardRuntimeTools(final Map params) {
        return STANDARD_RUNTIME_TOOLS.collect { tool ->
            [name: tool.name, path: runtimeToolPath(params, tool.binary), used: 'yes']
        }
    }

    static List<Map> configuredRuntimeTools(final Map params) {
        return standardRuntimeTools(params)
    }

    static Map runtimeContext(final Map params) {
        return [
            runtime_env_prefix: runtimeEnvPrefix(params),
            runtime_bin_dir   : runtimeBinDir(params),
            runtime_tmpdir    : runtimeTmpdir(params),
        ]
    }

    static void validateRuntimeContract(final Map params) {
        validateConfiguredDirectory('runtime env prefix', runtimeEnvPrefix(params))
        validateConfiguredDirectory('runtime bin dir', runtimeBinDir(params))
        validateConfiguredWritableDirectory('runtime tmpdir', runtimeTmpdir(params), true)

        standardRuntimeTools(params).each { tool ->
            validateConfiguredExecutable("runtime ${tool.name}", tool.path as String)
        }
    }

    static void validateConfiguredWritableDirectory(
        final String label,
        final String rawPath,
        final boolean createIfMissing = false
    ) {
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            throw new IllegalStateException("Missing configured writable directory path for ${label}")
        }

        final File directory = new File(path)
        if( !directory.exists() && createIfMissing ) {
            if( !directory.mkdirs() && !directory.exists() ) {
                throw new IllegalStateException(
                    "Configured writable directory for ${label} does not exist and could not be created: ${directory}"
                )
            }
        }

        if( !directory.exists() || !directory.isDirectory() ) {
            throw new IllegalStateException(
                "Configured writable directory for ${label} is missing or not a directory: ${directory}"
            )
        }
        if( !directory.canWrite() ) {
            throw new IllegalStateException(
                "Configured writable directory for ${label} is not writable: ${directory}"
            )
        }
    }

    static String shellExports(final Map params) {
        final String envPrefix = runtimeEnvPrefix(params)
        final String binDir = runtimeBinDir(params)
        final String tmpdir = runtimeTmpdir(params)
        final Map<String, String> exports = [
            RUNTIME_ENV_PREFIX     : envPrefix,
            RUNTIME_BIN_DIR        : binDir,
            TMPDIR                 : tmpdir,
            PYTHON3_BIN            : "${binDir}/python3",
            SAMTOOLS_BIN           : "${binDir}/samtools",
            BWA_MEM2_BIN           : "${binDir}/bwa-mem2",
            BAMCOVERAGE_BIN        : "${binDir}/bamCoverage",
            GATK_BIN               : "${binDir}/gatk",
        ]

        return exports.collect { key, value ->
            "export ${key}=${shellQuote(value)}"
        }.join('\n') + '\nmkdir -p "$TMPDIR"'
    }

    private static String shellQuote(final Object value) {
        final String text = (value ?: '').toString()
        return "'" + text.replace("'", "'\"'\"'") + "'"
    }

    static void writeRuntimeContract(
        final String rawOutdir,
        final List<Map> configuredTools,
        final Map runtimeContext = [:]
    ) {
        final File pipelineInfoDir = new File((rawOutdir ?: 'results').toString(), 'pipeline_info')
        if( !pipelineInfoDir.exists() ) {
            pipelineInfoDir.mkdirs()
        }

        final File reportFile = new File(pipelineInfoDir, 'runtime_contract.tsv')
        final StringBuilder builder = new StringBuilder()
        builder.append("tool\tconfigured_path\texists\tcurrently_used\n")

        configuredTools.each { tool ->
            final String path = (tool.path ?: '').toString()
            final boolean exists = path ? new File(path).exists() : false
            builder.append("${tool.name}\t${path}\t${exists}\t${tool.used}\n")
        }

        builder.append("\n[runtime_environment]\n")
        runtimeContext.each { key, value ->
            builder.append("${key}\t${(value ?: '').toString()}\n")
        }

        reportFile.text = builder.toString()
    }
}
