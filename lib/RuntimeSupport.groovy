class RuntimeSupport {

    private static final List<Map> STANDARD_RUNTIME_TOOLS = [
        [name: 'python3', binary: 'python3'],
        [name: 'trim_galore', binary: 'trim_galore'],
        [name: 'STAR', binary: 'STAR'],
        [name: 'samtools', binary: 'samtools'],
        [name: 'bedGraphToBigWig', binary: 'bedGraphToBigWig'],
        [name: 'bwa-mem2', binary: 'bwa-mem2'],
        [name: 'bamCoverage', binary: 'bamCoverage'],
        [name: 'gatk', binary: 'gatk'],
        [name: 'codon', binary: 'codon'],
        [name: 'pigz', binary: 'pigz'],
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

    static String runtimeCodonHome(final Map params) {
        return runtimeEnvPrefix(params)
    }

    static String resolveProjectPath(final String rawProjectDir, final Object rawPath) {
        final String projectDir = rawProjectDir?.toString()?.trim()
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            return path
        }

        final File candidate = new File(path)
        if( candidate.isAbsolute() || !projectDir ) {
            return candidate.path
        }

        return new File(projectDir, path).canonicalPath
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
            codon_home        : runtimeCodonHome(params),
            runtime_tmpdir    : runtimeTmpdir(params),
        ]
    }

    static void validateRuntimeContract(final Map params) {
        validateConfiguredDirectory('runtime env prefix', runtimeEnvPrefix(params))
        validateConfiguredDirectory('runtime bin dir', runtimeBinDir(params))
        validateConfiguredDirectory('codon home', runtimeCodonHome(params))
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
            TRIM_GALORE_BIN        : "${binDir}/trim_galore",
            STAR_BIN               : "${binDir}/STAR",
            SAMTOOLS_BIN           : "${binDir}/samtools",
            BEDGRAPH_TO_BIGWIG_BIN : "${binDir}/bedGraphToBigWig",
            BWA_MEM2_BIN           : "${binDir}/bwa-mem2",
            BAMCOVERAGE_BIN        : "${binDir}/bamCoverage",
            GATK_BIN               : "${binDir}/gatk",
            CODON_BIN              : "${binDir}/codon",
            PIGZ_BIN               : "${binDir}/pigz",
            CODON_HOME             : envPrefix,
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
        final String codonPreflightOutput,
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

        builder.append("\n[host_codon_seq_preflight]\n")
        builder.append(codonPreflightOutput ?: '')
        builder.append('\n')

        reportFile.text = builder.toString()
    }
}
