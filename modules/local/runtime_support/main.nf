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
        SAMTOOLS_BIN           : "${binDir}/samtools",
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
