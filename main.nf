#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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

def configuredRuntimeTools() {
    return [
        [name: 'python3', path: (params.runtime_python ?: '').toString(), used: 'yes'],
        [name: 'trim_galore', path: (params.runtime_trim_galore ?: '').toString(), used: 'yes'],
        [name: 'STAR', path: (params.runtime_star ?: '').toString(), used: 'yes'],
        [name: 'samtools', path: (params.runtime_samtools ?: '').toString(), used: 'yes'],
        [name: 'bedGraphToBigWig', path: (params.runtime_bedgraph_to_bigwig ?: '').toString(), used: 'yes'],
        [name: 'bwa-mem2', path: (params.runtime_bwa_mem2 ?: '').toString(), used: 'yes'],
        [name: 'bamCoverage', path: (params.runtime_bam_coverage ?: '').toString(), used: 'future'],
        [name: 'gatk', path: "${params.gatk_root ?: ''}/gatk".toString(), used: 'yes'],
    ]
}

def requireExecutable(final String label, final String rawPath) {
    final String path = rawPath?.toString()?.trim()
    if( !path ) {
        throw new IllegalStateException("Missing configured executable path for ${label}")
    }

    final File executable = new File(path)
    if( !executable.exists() || !executable.canExecute() ) {
        throw new IllegalStateException("Configured executable for ${label} is missing or not executable: ${executable}")
    }
}

def writeRuntimeContract(final String codonPreflightOutput) {
    final File pipelineInfoDir = new File((params.outdir ?: 'results').toString(), 'pipeline_info')
    if( !pipelineInfoDir.exists() ) {
        pipelineInfoDir.mkdirs()
    }

    final File reportFile = new File(pipelineInfoDir, 'runtime_contract.tsv')
    final StringBuilder builder = new StringBuilder()
    builder.append("tool\tconfigured_path\texists\tcurrently_used\n")

    configuredRuntimeTools().each { tool ->
        final String path = tool.path ?: ''
        final boolean exists = path ? new File(path).exists() : false
        builder.append("${tool.name}\t${path}\t${exists}\t${tool.used}\n")
    }

    builder.append("\n[host_codon_seq_preflight]\n")
    builder.append(codonPreflightOutput ?: '')
    builder.append('\n')

    reportFile.text = builder.toString()
}

requireExecutable('runtime python3', params.runtime_python as String)
final String codonPreflightOutput = enforcePinnedCodonSeq()
writeRuntimeContract(codonPreflightOutput)

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ()
}
