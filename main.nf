#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def enforcePinnedCodonSeq() {
    final File preflight = new File(projectDir.toString(), 'bin/check_codon_seq_host.sh')

    if( !preflight.exists() ) {
        throw new IllegalStateException(
            "Global pinned Codon/Seq preflight failed. Missing required preflight script: ${preflight}"
        )
    }

    final Process process = new ProcessBuilder('bash', preflight.toString())
        .directory(projectDir.toFile())
        .redirectErrorStream(true)
        .start()

    final String output = process.inputStream.getText('UTF-8').trim()
    final int exitCode = process.waitFor()

    if( exitCode != 0 ) {
        final String detail = output ? "\n${output}" : ''
        throw new IllegalStateException(
            "Global pinned Codon/Seq preflight failed. Every pipeline run requires codon 0.16.3 and Seq 0.11.3.${detail}"
        )
    }
}

enforcePinnedCodonSeq()

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ()
}
