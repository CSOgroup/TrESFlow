#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import PipelineSupport

def enforcePinnedCodonSeq() {
    final File preflight = new File(projectDir.toString(), 'bin/check_codon_seq_host.sh')
    final String runtimeEnvPrefix = (params.runtime_env_prefix ?: '').toString().trim()

    if( !preflight.exists() ) {
        throw new IllegalStateException(
            "Global pinned Codon/Seq preflight failed. Missing required preflight script: ${preflight}"
        )
    }

    final ProcessBuilder processBuilder = new ProcessBuilder('bash', preflight.toString())
        .directory(projectDir.toFile())
        .redirectErrorStream(true)

    if( runtimeEnvPrefix ) {
        final File runtimeBin = new File(runtimeEnvPrefix, 'bin')
        if( runtimeBin.exists() ) {
            final Map<String, String> env = processBuilder.environment()
            final String existingPath = env.get('PATH') ?: ''
            env.put('PATH', "${runtimeBin.toString()}:${existingPath}".toString())
        }
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

PipelineSupport.validateConfiguredExecutable('runtime python3', params.runtime_python as String)
final String codonPreflightOutput = enforcePinnedCodonSeq()
PipelineSupport.writeRuntimeContract(
    (params.outdir ?: 'results').toString(),
    PipelineSupport.configuredRuntimeTools(params),
    codonPreflightOutput
)

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ()
}
