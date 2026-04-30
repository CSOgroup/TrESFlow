class WorkflowSupport {

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

    static void validateReferenceContract(final Map references, final Map modalities) {
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

        final File ligationWhitelist = new File(references.ligation_barcode_whitelist as String)
        if( !ligationWhitelist.exists() || !ligationWhitelist.isFile() ) {
            throw new IllegalArgumentException(
                "Required ligation barcode whitelist not found: ${ligationWhitelist}. " +
                "Expected <references.root>/ligation_barcode_whitelist.txt"
            )
        }

        if( (modalities.rna as boolean) ) {
            validateRnaReferences(
                references.rna_star_index_dir as String,
                references.rna_chrom_sizes as String,
                root
            )
        }

        if( (modalities.dna as boolean) ) {
            validateDnaReferences(
                references.dna_bwa_reference as String,
                references.dna_blacklist_bed as String,
                references.dna_effective_genome_size_file as String,
                root
            )
        }
    }

    static String requireBwaMem2Prefix(final String rawPrefix) {
        final String prefix = rawPrefix?.toString()?.trim()
        if( !prefix ) {
            throw new IllegalArgumentException(
                "Missing required DNA bwa-mem2 prefix resolved from samplesheet field references.root"
            )
        }

        final File prefixFile = new File(prefix)
        if( !prefixFile.exists() ) {
            throw new IllegalArgumentException(
                "DNA bwa-mem2 prefix FASTA not found: ${prefixFile}. " +
                "Expected <references.root>/dna/human/bwa/hg38.fa"
            )
        }

        final List<String> suffixes = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
        suffixes.each { suffix ->
            final File sidecar = new File("${prefix}${suffix}")
            if( !sidecar.exists() ) {
                throw new IllegalArgumentException(
                    "DNA bwa-mem2 index sidecar not found for prefix '${prefix}': ${sidecar}. " +
                    "Controlled by samplesheet field references.root"
                )
            }
        }

        return prefix
    }

    static void validateRnaReferences(
        final String rawStarIndexDir,
        final String rawChromSizes,
        final String rawReferenceRoot = ''
    ) {
        final String starIndexDir = rawStarIndexDir?.toString()?.trim()
        final String chromSizes = rawChromSizes?.toString()?.trim()

        if( !starIndexDir ) {
            throw new IllegalArgumentException(
                "Missing required RNA STAR index directory resolved from samplesheet field references.root"
            )
        }

        final File starDir = new File(starIndexDir)
        if( !starDir.exists() || !starDir.isDirectory() ) {
            throw new IllegalArgumentException(
                "RNA STAR index directory not found: ${starDir}. " +
                "Expected <references.root>/rna/human/star"
            )
        }

        if( !chromSizes ) {
            throw new IllegalArgumentException(
                "Missing required RNA chromosome sizes file resolved from samplesheet field references.root"
            )
        }

        final File chromSizesFile = new File(chromSizes)
        if( !chromSizesFile.exists() || !chromSizesFile.isFile() ) {
            throw new IllegalArgumentException(
                "RNA chromosome sizes file not found: ${chromSizesFile}. " +
                "Expected <references.root>/rna/human/chrom.sizes"
            )
        }
    }

    static void validateRnaAlignment(final String rawRefBaseDir, final String rawSpecies) {
        throw new IllegalArgumentException(
            "Legacy RNA reference parameters are no longer supported. Use samplesheet field references.root"
        )
    }

    static long validateDnaReferences(
        final String rawBwaReference,
        final String rawBlacklistBed,
        final String rawEffectiveGenomeSizeFile,
        final String rawReferenceRoot = ''
    ) {
        requireBwaMem2Prefix(rawBwaReference)

        final String blacklistBed = rawBlacklistBed?.toString()?.trim()
        final String effectiveGenomeSizeFile = rawEffectiveGenomeSizeFile?.toString()?.trim()

        if( !blacklistBed ) {
            throw new IllegalArgumentException(
                "Missing required DNA blacklist BED resolved from samplesheet field references.root"
            )
        }

        final File blacklistFile = new File(blacklistBed)
        if( !blacklistFile.exists() ) {
            throw new IllegalArgumentException(
                "DNA blacklist BED not found: ${blacklistFile}. Expected <references.root>/dna/human/blacklist.bed"
            )
        }

        if( !effectiveGenomeSizeFile ) {
            throw new IllegalArgumentException(
                "Missing required DNA effective genome size file resolved from samplesheet field references.root"
            )
        }

        final File effSizeFile = new File(effectiveGenomeSizeFile)
        if( !effSizeFile.exists() || !effSizeFile.isFile() ) {
            throw new IllegalArgumentException(
                "DNA effective genome size file not found: ${effSizeFile}. " +
                "Expected <references.root>/dna/human/effective_genome_size.txt"
            )
        }

        final String effSizeRaw = effSizeFile.text.trim()
        long effSize = 0L
        try {
            effSize = effSizeRaw as long
        }
        catch( Exception ignored ) {
            throw new IllegalArgumentException(
                "Invalid DNA effective genome size in ${effSizeFile}: '${effSizeRaw}'. Value must be an integer > 0"
            )
        }

        if( effSize < 1 ) {
            throw new IllegalArgumentException(
                "Invalid DNA effective genome size in ${effSizeFile}: '${effSizeRaw}'. Value must be an integer > 0"
            )
        }

        return effSize
    }

    static void validateDnaAlignment(
        final String rawBwaReference,
        final String rawBlacklistBed,
        final String rawEffectiveGenomeSize
    ) {
        throw new IllegalArgumentException(
            "Legacy DNA reference parameters are no longer supported. Use samplesheet field references.root"
        )
    }
}
