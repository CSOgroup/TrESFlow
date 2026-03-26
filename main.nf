#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import RuntimeSupport

// The runtime contract is enforced once up front so every downstream task sees
// the same validated toolchain and the same pinned Codon/Seq preflight result.
def runCodonSeqPreflight() {
    final File preflight = new File(projectDir.toString(), 'bin/check_codon_seq_host.sh')
    final String codonBin = RuntimeSupport.runtimeToolPath(params, 'codon')
    final String codonHome = RuntimeSupport.runtimeCodonHome(params)

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

RuntimeSupport.validateRuntimeContract(params)
final String codonPreflightOutput = runCodonSeqPreflight()
RuntimeSupport.writeRuntimeContract(
    (params.outdir ?: 'results').toString(),
    RuntimeSupport.configuredRuntimeTools(params),
    codonPreflightOutput,
    RuntimeSupport.runtimeContext(params)
)

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ()
}
