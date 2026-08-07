#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CANONICAL_BAM_FILTER_TEST {
    tag 'ucsc-and-ensembl'
    label 'process_single'

    input:
    val repoRoot

    output:
    path 'canonical_bam_filter.passed'

    script:
    """
    TRESFLOW_REPO_ROOT="${repoRoot}" \
      bash "${repoRoot}/tests/test_canonical_bam_filter.sh"
    touch canonical_bam_filter.passed
    """
}

workflow {
    def repoRoot = params.repo_root ?: new java.io.File(projectDir.toString(), '..').canonicalPath
    CANONICAL_BAM_FILTER_TEST(repoRoot)
}
