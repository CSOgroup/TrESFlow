#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import PipelineSupport

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
    'runtime python3': params.runtime_python,
    'runtime trim_galore': params.runtime_trim_galore,
    'runtime STAR': params.runtime_star,
    'runtime samtools': params.runtime_samtools,
    'runtime bedGraphToBigWig': params.runtime_bedgraph_to_bigwig,
    'runtime bwa-mem2': params.runtime_bwa_mem2,
    'runtime bamCoverage': params.runtime_bam_coverage,
    'runtime gatk': params.runtime_gatk,
    'runtime codon': params.runtime_codon,
].each { label, path ->
    PipelineSupport.validateConfiguredExecutable(label, path as String)
}
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
