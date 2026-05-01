class WorkflowSupport {

    private static final List<String> REQUIRED_STAR_INDEX_FILES = [
        'Genome',
        'SA',
        'SAindex',
        'chrName.txt',
        'chrLength.txt',
        'chrStart.txt',
        'chrNameLength.txt',
        'genomeParameters.txt',
    ]

    private static final List<String> REQUIRED_BWA_MEM2_SIDECARS = [
        '.0123',
        '.amb',
        '.ann',
        '.bwt.2bit.64',
        '.pac',
    ]

    static List asPathList(final Object value) {
        if( value instanceof List ) {
            return value
        }
        return [value]
    }

    static List<Map> pairRnaSplitFastqs(final String sampleId, final Object splitR1s, final Object splitR2s) {
        final Map<String, Object> r1ByGroup = asPathList(splitR1s).collectEntries { path ->
            final String name = path.getName()
            final String group = name
                .replaceFirst("^${sampleId}_", '')
                .replaceFirst('_R1\\.(?:fastq|fq)\\.gz$', '')
            [(group): path]
        }
        final Map<String, Object> r2ByGroup = asPathList(splitR2s).collectEntries { path ->
            final String name = path.getName()
            final String group = name
                .replaceFirst("^${sampleId}_", '')
                .replaceFirst('_R2\\.(?:fastq|fq)\\.gz$', '')
            [(group): path]
        }

        final List<String> groups = (r1ByGroup.keySet() + r2ByGroup.keySet()).unique().sort()
        return groups.collect { group ->
            if( !r1ByGroup.containsKey(group) || !r2ByGroup.containsKey(group) ) {
                throw new IllegalStateException(
                    "Missing split FASTQ mate for sample '${sampleId}' group '${group}'"
                )
            }

            [
                splitName: "${sampleId}_${group}",
                r1       : r1ByGroup[group],
                r2       : r2ByGroup[group],
            ]
        }
    }

    static List<Map> pairDnaSplitFastqs(final Object splitR1s, final Object splitR2s) {
        final Map<String, Object> r1BySplit = asPathList(splitR1s).collectEntries { path ->
            [(path.getName().replaceFirst('_R1\\.(?:fastq|fq)\\.gz$', '')): path]
        }
        final Map<String, Object> r2BySplit = asPathList(splitR2s).collectEntries { path ->
            [(path.getName().replaceFirst('_R2\\.(?:fastq|fq)\\.gz$', '')): path]
        }

        final List<String> splitNames = (r1BySplit.keySet() + r2BySplit.keySet()).unique().sort()
        return splitNames.collect { splitName ->
            if( !r1BySplit.containsKey(splitName) || !r2BySplit.containsKey(splitName) ) {
                throw new IllegalStateException("Missing split FASTQ mate for DNA split '${splitName}'")
            }

            [
                splitName: splitName,
                r1       : r1BySplit[splitName],
                r2       : r2BySplit[splitName],
            ]
        }
    }

    static List<Map> collectDnaRgHeaders(final Object rgHeaders) {
        return asPathList(rgHeaders).collect { rgHeader ->
            [
                splitName: rgHeader.getName()
                    .replaceFirst('^SAM_RG_Header_', '')
                    .replaceFirst('\\.tsv$', ''),
                rgHeader : rgHeader,
            ]
        }
    }

    static Map parseDnaSplitName(final String sampleId, final String splitName) {
        final String suffix = splitName.replaceFirst("^${sampleId}_", '')
        final List<String> tokens = suffix.tokenize('_')
        if( tokens.size() < 2 ) {
            throw new IllegalStateException(
                "Unable to derive DNA group and modality from split output '${splitName}'"
            )
        }

        final String group = tokens[0]
        return [
            group     : group,
            modality  : tokens[1..-1].join('_'),
            sampleGroup: "${sampleId}_${group}",
        ]
    }

    static void validateReferenceContract(
        final Map references,
        final Map modalities,
        final List<Map> sampleRows = []
    ) {
        requireDirectory(
            references.root,
            'references.root',
            "Reference root not found or not a directory"
        )
        requireValue(references.species, 'references.species')
        requireFile(
            references.ligation_barcode_whitelist,
            'references.ligation_barcode_whitelist',
            "Required ligation barcode whitelist not found"
        )

        if( (modalities.rna as boolean) ) {
            validateRnaReferences(
                references.rna_ref_dir as String,
                references.rna_chrom_sizes as String
            )
        }

        if( (modalities.dna as boolean) ) {
            final String bwaPrefix = validateDnaReferences(
                references.dna_ref_dir as String,
                references.dna_blacklist_bed as String,
                references.dna_effective_genome_size as String,
                references.dna_chrom_sizes as String
            )
            references.dna_bwa_reference = bwaPrefix
            sampleRows.findAll { row -> row.modality == 'dna' }.each { row ->
                row.dna_bwa_reference = bwaPrefix
            }
        }
    }

    static String inferBwaMem2Prefix(final String rawDnaRefDir) {
        final File dir = requireDirectory(
            rawDnaRefDir,
            'references.dna_ref_dir',
            "DNA bwa-mem2 index directory not found or not a directory"
        )

        final Map<String, Set<String>> suffixesByPrefix = [:].withDefault { new LinkedHashSet<String>() }
        dir.listFiles()?.findAll { file -> file.isFile() }?.each { file ->
            REQUIRED_BWA_MEM2_SIDECARS.each { suffix ->
                if( file.name.endsWith(suffix) ) {
                    final String prefixName = file.name.substring(0, file.name.length() - suffix.length())
                    suffixesByPrefix[new File(dir, prefixName).canonicalPath].add(suffix)
                }
            }
        }

        final List<String> completePrefixes = suffixesByPrefix.findAll { prefix, suffixes ->
            REQUIRED_BWA_MEM2_SIDECARS.every { suffix -> suffixes.contains(suffix) }
        }.keySet().sort()

        if( completePrefixes.size() == 1 ) {
            return completePrefixes[0]
        }

        if( completePrefixes.size() > 1 ) {
            throw new IllegalArgumentException(
                "Multiple complete bwa-mem2 sidecar sets found in references.dna_ref_dir '${dir}': " +
                completePrefixes.join(', ') + ". Keep exactly one complete index prefix in this directory."
            )
        }

        if( suffixesByPrefix.size() == 1 ) {
            final String prefix = suffixesByPrefix.keySet().first()
            final List<String> missing = REQUIRED_BWA_MEM2_SIDECARS.findAll { suffix ->
                !suffixesByPrefix[prefix].contains(suffix)
            }
            throw new IllegalArgumentException(
                "No complete bwa-mem2 sidecar set found in references.dna_ref_dir '${dir}'. " +
                "Detected prefix '${prefix}' is missing: " + missing.collect { suffix -> "${prefix}${suffix}" }.join(', ')
            )
        }

        throw new IllegalArgumentException(
            "No bwa-mem2 sidecar files found in references.dna_ref_dir '${dir}'. " +
            "Expected exactly one complete set: " + REQUIRED_BWA_MEM2_SIDECARS.collect { suffix -> "<prefix>${suffix}" }.join(', ')
        )
    }

    static String requireBwaMem2Prefix(final String rawPrefix) {
        final String prefix = rawPrefix?.toString()?.trim()
        if( !prefix ) {
            throw new IllegalArgumentException(
                "Missing required DNA bwa-mem2 prefix"
            )
        }

        REQUIRED_BWA_MEM2_SIDECARS.each { suffix ->
            final File sidecar = new File("${prefix}${suffix}")
            if( !sidecar.exists() ) {
                throw new IllegalArgumentException(
                    "DNA bwa-mem2 index sidecar not found for inferred prefix '${prefix}': ${sidecar}"
                )
            }
        }

        return prefix
    }

    static void validateRnaReferences(
        final String rawStarIndexDir,
        final String rawChromSizes
    ) {
        final File starDir = requireDirectory(
            rawStarIndexDir,
            'references.rna_ref_dir',
            "RNA STAR index directory not found"
        )

        REQUIRED_STAR_INDEX_FILES.each { fileName ->
            requireFile(
                new File(starDir, fileName).path,
                'references.rna_ref_dir',
                "RNA STAR index file missing for references.rna_ref_dir"
            )
        }

        requireFile(
            rawChromSizes,
            'references.rna_ref_dir',
            "RNA chromosome sizes file for coverage not found. Expected chrNameLength.txt inside the STAR index directory"
        )
    }

    static String validateDnaReferences(
        final String rawDnaRefDir,
        final String rawBlacklistBed,
        final String rawEffectiveGenomeSize,
        final String rawDnaChromSizes = ''
    ) {
        final String bwaPrefix = inferBwaMem2Prefix(rawDnaRefDir)
        requireBwaMem2Prefix(bwaPrefix)

        final String effSizeRaw = rawEffectiveGenomeSize?.toString()?.trim()
        final String dnaChromSizes = rawDnaChromSizes?.toString()?.trim()

        requireFile(rawBlacklistBed, 'references.dna_blacklist_bed', "DNA blacklist BED not found")

        if( dnaChromSizes ) {
            requireFile(dnaChromSizes, 'references.dna_chrom_sizes', "DNA chromosome sizes file not found")
        }

        if( !effSizeRaw ) {
            throw new IllegalArgumentException(
                "Missing required DNA effective genome size: references.dna_effective_genome_size. " +
                "This is used by BAM_COVERAGE_DNA for bamCoverage --effectiveGenomeSize."
            )
        }

        long effSize = 0L
        try {
            effSize = effSizeRaw as long
        }
        catch( Exception ignored ) {
            throw new IllegalArgumentException(
                "Invalid references.dna_effective_genome_size '${effSizeRaw}'. Value must be an integer > 0"
            )
        }

        if( effSize < 1 ) {
            throw new IllegalArgumentException(
                "Invalid references.dna_effective_genome_size '${effSizeRaw}'. Value must be an integer > 0"
            )
        }

        return bwaPrefix
    }

    private static String requireValue(final Object rawValue, final String fieldName) {
        final String value = rawValue?.toString()?.trim()
        if( !value ) {
            throw new IllegalArgumentException("Missing required samplesheet field: ${fieldName}")
        }
        return value
    }

    private static File requireDirectory(final Object rawPath, final String fieldName, final String message) {
        final File directory = new File(requireValue(rawPath, fieldName))
        if( !directory.exists() || !directory.isDirectory() ) {
            throw new IllegalArgumentException(
                "${message}: ${directory}. Controlled by samplesheet field ${fieldName}"
            )
        }
        return directory
    }

    private static File requireFile(final Object rawPath, final String fieldName, final String message) {
        final File file = new File(requireValue(rawPath, fieldName))
        if( !file.exists() || !file.isFile() ) {
            throw new IllegalArgumentException(
                "${message}: ${file}. Controlled by samplesheet field ${fieldName}"
            )
        }
        return file
    }
}
