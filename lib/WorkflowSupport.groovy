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
                .replaceFirst('_R1\\.fq\\.gz$', '')
            [(group): path]
        }
        final Map<String, Object> r2ByGroup = asPathList(splitR2s).collectEntries { path ->
            final String name = path.getName()
            final String group = name
                .replaceFirst("^${sampleId}_", '')
                .replaceFirst('_R2\\.fq\\.gz$', '')
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
            [(path.getName().replaceFirst('_R1\\.fq\\.gz$', '')): path]
        }
        final Map<String, Object> r2BySplit = asPathList(splitR2s).collectEntries { path ->
            [(path.getName().replaceFirst('_R2\\.fq\\.gz$', '')): path]
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
        final String root = references.root?.toString()?.trim()
        if( !root ) {
            throw new IllegalArgumentException("Missing required samplesheet field: references.root")
        }

        final File rootDir = new File(root)
        if( !rootDir.exists() || !rootDir.isDirectory() ) {
            throw new IllegalArgumentException(
                "Reference root not found or not a directory: ${rootDir}. Controlled by samplesheet field references.root"
            )
        }

        final String species = references.species?.toString()?.trim()
        if( !species ) {
            throw new IllegalArgumentException("Missing required samplesheet field: references.species")
        }

        final File ligationWhitelist = new File(references.ligation_barcode_whitelist as String)
        if( !ligationWhitelist.exists() || !ligationWhitelist.isFile() ) {
            throw new IllegalArgumentException(
                "Required ligation barcode whitelist not found: ${ligationWhitelist}. " +
                "Controlled by samplesheet field references.ligation_barcode_whitelist"
            )
        }

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
        final String dnaRefDir = rawDnaRefDir?.toString()?.trim()
        if( !dnaRefDir ) {
            throw new IllegalArgumentException(
                "Missing required DNA bwa-mem2 index directory: references.dna_ref_dir"
            )
        }

        final File dir = new File(dnaRefDir)
        if( !dir.exists() || !dir.isDirectory() ) {
            throw new IllegalArgumentException(
                "DNA bwa-mem2 index directory not found or not a directory: ${dir}. " +
                "Controlled by samplesheet field references.dna_ref_dir"
            )
        }

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
        final String starIndexDir = rawStarIndexDir?.toString()?.trim()
        final String chromSizes = rawChromSizes?.toString()?.trim()

        if( !starIndexDir ) {
            throw new IllegalArgumentException(
                "Missing required RNA STAR index directory: references.rna_ref_dir"
            )
        }

        final File starDir = new File(starIndexDir)
        if( !starDir.exists() || !starDir.isDirectory() ) {
            throw new IllegalArgumentException(
                "RNA STAR index directory not found for references.rna_ref_dir: ${starDir}"
            )
        }

        REQUIRED_STAR_INDEX_FILES.each { fileName ->
            final File requiredFile = new File(starDir, fileName)
            if( !requiredFile.exists() || !requiredFile.isFile() ) {
                throw new IllegalArgumentException(
                    "RNA STAR index file missing for references.rna_ref_dir: ${requiredFile}"
                )
            }
        }

        if( !chromSizes ) {
            throw new IllegalArgumentException(
                "Missing RNA chromosome sizes file derived from references.rna_ref_dir"
            )
        }
        final File chromSizesFile = new File(chromSizes)
        if( !chromSizesFile.exists() || !chromSizesFile.isFile() ) {
            throw new IllegalArgumentException(
                "RNA chromosome sizes file for coverage not found for references.rna_ref_dir: ${chromSizesFile}. " +
                "Expected chrNameLength.txt inside the STAR index directory"
            )
        }
    }

    static void validateRnaAlignment(final String rawRefBaseDir, final String rawSpecies) {
        throw new IllegalArgumentException(
            "Legacy RNA reference parameters are no longer supported. Use samplesheet field references.rna_ref_dir"
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

        final String blacklistBed = rawBlacklistBed?.toString()?.trim()
        final String effSizeRaw = rawEffectiveGenomeSize?.toString()?.trim()
        final String dnaChromSizes = rawDnaChromSizes?.toString()?.trim()

        if( !blacklistBed ) {
            throw new IllegalArgumentException(
                "Missing required DNA blacklist BED: references.dna_blacklist_bed"
            )
        }

        final File blacklistFile = new File(blacklistBed)
        if( !blacklistFile.exists() ) {
            throw new IllegalArgumentException(
                "DNA blacklist BED not found for references.dna_blacklist_bed: ${blacklistFile}"
            )
        }

        if( dnaChromSizes ) {
            final File chromSizesFile = new File(dnaChromSizes)
            if( !chromSizesFile.exists() || !chromSizesFile.isFile() ) {
                throw new IllegalArgumentException(
                    "DNA chromosome sizes file not found for references.dna_chrom_sizes: ${chromSizesFile}"
                )
            }
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

    static void validateDnaAlignment(
        final String rawBwaReference,
        final String rawBlacklistBed,
        final String rawEffectiveGenomeSize
    ) {
        throw new IllegalArgumentException(
            "Legacy DNA reference parameters are no longer supported. Use samplesheet fields references.dna_ref_dir, references.dna_blacklist_bed, and references.dna_effective_genome_size"
        )
    }
}
