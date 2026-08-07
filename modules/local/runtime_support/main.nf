/*
 * Strict-syntax-compatible runtime environment rendering shared by local
 * processes. Full host-side validation and reporting remain in lib/RuntimeSupport.groovy.
 */

def runtimeShellExports(runtimeParams) {
    def envPrefix = (runtimeParams.runtime_env_prefix ?: '').toString().trim()
    def binDir = envPrefix ? "${envPrefix}/bin" : ''
    def tmpdir = (runtimeParams.runtime_tmpdir ?: '').toString().trim()
    def exports = [
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

    exports.collect { key, value ->
        def text = (value ?: '').toString()
        def quoted = "'" + text.replace("'", "'\"'\"'") + "'"
        "export ${key}=${quoted}"
    }.join('\n') + '\nmkdir -p "$TMPDIR"'
}

def runtimeOutdir() {
    def value = java.lang.System.getProperty('tresflow.resolvedOutdir')
    if( !value ) {
        throw new IllegalStateException('TrESFlow resolved output directory is not initialized')
    }
    return value
}

def runtimeCoreScriptsDir() {
    def value = java.lang.System.getProperty('tresflow.resolvedCoreScriptsDir')
    if( !value ) {
        throw new IllegalStateException('TrESFlow resolved core scripts directory is not initialized')
    }
    return value
}
