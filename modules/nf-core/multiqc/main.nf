process MULTIQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/biocontainers/multiqc@sha256:dfd9fde2c48b896b884e79a71ddc16c72c97a0ee5c5c8e45aaba50f55d07d263'
        : 'quay.io/biocontainers/multiqc@sha256:dfd9fde2c48b896b884e79a71ddc16c72c97a0ee5c5c8e45aaba50f55d07d263'}"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*"), path(multiqc_config, stageAs: "?/*"), path(multiqc_logo), path(replace_names), path(sample_names)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), emit: plots, optional: true
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? multiqc_config instanceof List ? "--config ${multiqc_config.join(' --config ')}" : "--config ${multiqc_config}" : ""
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    export TMPDIR="\$PWD/.tmp"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export MPLCONFIGDIR="\$PWD/.mplconfig"
    export NUMBA_CACHE_DIR="\$PWD/.numba"
    mkdir -p "\$TMPDIR" "\$XDG_CACHE_HOME" "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

    multiqc \\
        --force \\
        --ignore '.tmp/**' \\
        --ignore '.cache/**' \\
        --ignore '.mplconfig/**' \\
        --ignore '.numba/**' \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}
