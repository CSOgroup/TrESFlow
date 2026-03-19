#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRESEQ } from './workflows/treseq'

workflow {
    TRESEQ()
}
