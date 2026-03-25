class WorkflowSupport {

    static List asPathList(final Object value) {
        if( value instanceof List ) {
            return value
        }
        return [value]
    }

    static String requireBwaMem2Prefix(final String rawPrefix) {
        final String prefix = rawPrefix?.toString()?.trim()
        if( !prefix ) {
            throw new IllegalArgumentException(
                "Missing required DNA alignment resource: resources.dna_bwa_reference or --dna_bwa_reference"
            )
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
                "Missing required RNA alignment resource: resources.rna_align_species or --rna_align_species (human|mouse)"
            )
        }
        if( !(species in ['human', 'mouse']) ) {
            throw new IllegalArgumentException(
                "Invalid RNA align species '${species}'. Supported values: human, mouse"
            )
        }
        if( !refBaseDir ) {
            throw new IllegalArgumentException(
                "Missing required RNA alignment resource: resources.rna_ref_base_dir or --rna_ref_base_dir"
            )
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
            throw new IllegalArgumentException(
                "Missing required DNA alignment resource: resources.dna_blacklist_bed or --dna_blacklist_bed"
            )
        }

        final File blacklistFile = new File(blacklistBed)
        if( !blacklistFile.exists() ) {
            throw new IllegalArgumentException("DNA blacklist BED not found: ${blacklistFile}")
        }

        if( !effSizeRaw ) {
            throw new IllegalArgumentException(
                "Missing required DNA alignment resource: resources.dna_effective_genome_size or --dna_effective_genome_size"
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
}
