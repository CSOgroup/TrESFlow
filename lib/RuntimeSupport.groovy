class RuntimeSupport {

    static void validateConfiguredExecutable(final String label, final String rawPath) {
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            throw new IllegalStateException("Missing configured executable path for ${label}")
        }

        final File executable = new File(path)
        if( !executable.exists() || !executable.canExecute() ) {
            throw new IllegalStateException(
                "Configured executable for ${label} is missing or not executable: ${executable}"
            )
        }
    }

    static void validateConfiguredDirectory(final String label, final String rawPath) {
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            throw new IllegalStateException("Missing configured directory path for ${label}")
        }

        final File directory = new File(path)
        if( !directory.exists() || !directory.isDirectory() ) {
            throw new IllegalStateException(
                "Configured directory for ${label} is missing or not a directory: ${directory}"
            )
        }
    }

    static List<Map> configuredRuntimeTools(final Map params) {
        return [
            [name: 'python3', path: (params.runtime_python ?: '').toString(), used: 'yes'],
            [name: 'trim_galore', path: (params.runtime_trim_galore ?: '').toString(), used: 'yes'],
            [name: 'STAR', path: (params.runtime_star ?: '').toString(), used: 'yes'],
            [name: 'samtools', path: (params.runtime_samtools ?: '').toString(), used: 'yes'],
            [name: 'bedGraphToBigWig', path: (params.runtime_bedgraph_to_bigwig ?: '').toString(), used: 'yes'],
            [name: 'bwa-mem2', path: (params.runtime_bwa_mem2 ?: '').toString(), used: 'yes'],
            [name: 'bamCoverage', path: (params.runtime_bam_coverage ?: '').toString(), used: 'yes'],
            [name: 'gatk', path: (params.runtime_gatk ?: '').toString(), used: 'yes'],
            [name: 'codon', path: (params.runtime_codon ?: '').toString(), used: 'yes'],
        ]
    }

    static void writeRuntimeContract(
        final String rawOutdir,
        final List<Map> configuredTools,
        final String codonPreflightOutput
    ) {
        final File pipelineInfoDir = new File((rawOutdir ?: 'results').toString(), 'pipeline_info')
        if( !pipelineInfoDir.exists() ) {
            pipelineInfoDir.mkdirs()
        }

        final File reportFile = new File(pipelineInfoDir, 'runtime_contract.tsv')
        final StringBuilder builder = new StringBuilder()
        builder.append("tool\tconfigured_path\texists\tcurrently_used\n")

        configuredTools.each { tool ->
            final String path = (tool.path ?: '').toString()
            final boolean exists = path ? new File(path).exists() : false
            builder.append("${tool.name}\t${path}\t${exists}\t${tool.used}\n")
        }

        builder.append("\n[host_codon_seq_preflight]\n")
        builder.append(codonPreflightOutput ?: '')
        builder.append('\n')

        reportFile.text = builder.toString()
    }
}
