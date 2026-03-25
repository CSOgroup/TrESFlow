#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import RuntimeSupport

def enforcePinnedCodonSeq() {
    final File preflight = new File(projectDir.toString(), 'bin/check_codon_seq_host.sh')
    final String codonBin = (params.runtime_codon ?: '').toString().trim()
    final String codonHome = (params.codon_home ?: '').toString().trim()

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

[
    ['runtime environment prefix', RuntimeSupport.runtimeEnvPrefix(params)],
    ['runtime bin dir', RuntimeSupport.runtimeBinDir(params)],
].each { label, path ->
    RuntimeSupport.validateConfiguredDirectory(label as String, path as String)
}

RuntimeSupport.standardRuntimeTools(params).each { tool ->
    RuntimeSupport.validateConfiguredExecutable("runtime ${tool.name}", tool.path as String)
}
RuntimeSupport.validateConfiguredExecutable('runtime codon', params.runtime_codon as String)
final String codonPreflightOutput = enforcePinnedCodonSeq()
RuntimeSupport.writeRuntimeContract(
    (params.outdir ?: 'results').toString(),
    RuntimeSupport.configuredRuntimeTools(params),
    codonPreflightOutput,
    [
        runtime_env_prefix: RuntimeSupport.runtimeEnvPrefix(params),
        runtime_bin_dir   : RuntimeSupport.runtimeBinDir(params),
    ]
)

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ()
}
