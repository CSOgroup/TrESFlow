class PipelineSupport {

    static List asPathList(final Object value) {
        if( value instanceof List ) {
            return value
        }
        return [value]
    }

    static String requireBwaMem2Prefix(final String rawPrefix) {
        final String prefix = rawPrefix?.toString()?.trim()
        if( !prefix ) {
            throw new IllegalArgumentException("Missing required parameter for DNA alignment: --dna_bwa_reference")
        }

        final List<String> suffixes = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
        suffixes.each { suffix ->
            final File sidecar = new File("${prefix}${suffix}")
            if( !sidecar.exists() ) {
                throw new IllegalArgumentException(
                    "DNA bwa-mem2 index sidecar not found for prefix '${prefix}': ${sidecar}"
                )
            }
        }

        return prefix
    }

    static void validateRnaAlignment(final String rawRefBaseDir, final String rawSpecies) {
        final String species = rawSpecies?.toString()?.trim()?.toLowerCase()
        final String refBaseDir = rawRefBaseDir?.toString()?.trim()

        if( !species ) {
            throw new IllegalArgumentException(
                "Missing required parameter for RNA alignment: --rna_align_species human|mouse"
            )
        }
        if( !(species in ['human', 'mouse']) ) {
            throw new IllegalArgumentException(
                "Invalid --rna_align_species '${species}'. Supported values: human, mouse"
            )
        }
        if( !refBaseDir ) {
            throw new IllegalArgumentException("Missing required parameter for RNA alignment: --rna_ref_base_dir")
        }

        final List<String> requiredRefPaths = species == 'human'
            ? ['GRCh38_TrES/star', 'hg38.chrom.sizes']
            : ['GRCm39_TrES/star', 'mm39.chrom.sizes']

        requiredRefPaths.each { relPath ->
            final File resolved = new File(refBaseDir, relPath)
            if( !resolved.exists() ) {
                throw new IllegalArgumentException(
                    "RNA alignment reference path not found for species '${species}': ${resolved}"
                )
            }
        }
    }

    static void validateDnaAlignment(
        final String rawBwaReference,
        final String rawBlacklistBed,
        final String rawEffectiveGenomeSize
    ) {
        requireBwaMem2Prefix(rawBwaReference)

        final String blacklistBed = rawBlacklistBed?.toString()?.trim()
        final String effSizeRaw = rawEffectiveGenomeSize?.toString()?.trim()

        if( !blacklistBed ) {
            throw new IllegalArgumentException("Missing required parameter for DNA alignment: --dna_blacklist_bed")
        }

        final File blacklistFile = new File(blacklistBed)
        if( !blacklistFile.exists() ) {
            throw new IllegalArgumentException("DNA blacklist BED not found: ${blacklistFile}")
        }

        if( !effSizeRaw ) {
            throw new IllegalArgumentException(
                "Missing required parameter for DNA alignment: --dna_effective_genome_size"
            )
        }

        long effSize = 0L
        try {
            effSize = effSizeRaw as long
        }
        catch( Exception ignored ) {
            throw new IllegalArgumentException(
                "Invalid --dna_effective_genome_size '${effSizeRaw}'. Value must be an integer > 0"
            )
        }

        if( effSize < 1 ) {
            throw new IllegalArgumentException(
                "Invalid --dna_effective_genome_size '${effSizeRaw}'. Value must be an integer > 0"
            )
        }
    }

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
