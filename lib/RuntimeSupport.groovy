class RuntimeSupport {

    private static final List<String> AUTOSOMES = (1..22).collect { it.toString() }
    private static final Set<String> UCSC_NAMES = (
        AUTOSOMES.collect { chromosome -> 'chr' + chromosome } + ['chrX', 'chrY', 'chrM']
    ) as Set<String>
    private static final Set<String> ENSEMBL_NAMES = (
        AUTOSOMES + ['X', 'Y', 'MT', 'M']
    ) as Set<String>

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

    static String resolvePath(final String rawBaseDir, final Object rawPath) {
        final String baseDir = rawBaseDir?.toString()?.trim()
        final String path = rawPath?.toString()?.trim()
        if( !path ) {
            return path
        }

        final File candidate = new File(path)
        if( candidate.isAbsolute() || !baseDir ) {
            return candidate.canonicalPath
        }

        return new File(baseDir, path).canonicalPath
    }

    static String resolveLaunchPath(final String rawLaunchDir, final Object rawPath) {
        return resolvePath(rawLaunchDir, rawPath)
    }

    static String resolveProjectPath(final String rawProjectDir, final Object rawPath) {
        return resolvePath(rawProjectDir, rawPath)
    }

    static String resolvePipelineReleaseVersion(
        final String rawProjectDir,
        final Object manifestVersion
    ) {
        final File projectDirectory = new File(rawProjectDir).canonicalFile
        final File resolver = new File(projectDirectory, 'bin/resolve_tresflow_release_version.sh')
        if( !resolver.exists() ) {
            throw new IllegalStateException(
                "Missing repository release-version resolver: ${resolver}"
            )
        }

        final Process process = new ProcessBuilder(
            'bash',
            resolver.canonicalPath,
            projectDirectory.canonicalPath,
            (manifestVersion ?: '').toString()
        )
            .directory(projectDirectory)
            .redirectErrorStream(true)
            .start()
        final String output = process.inputStream.getText('UTF-8').trim()
        final int exitCode = process.waitFor()
        if( exitCode != 0 || !output ) {
            throw new IllegalStateException(
                "Unable to resolve the TrESFlow release version: ${output ?: 'no output'}"
            )
        }
        return output
    }

    static Map writeCanonicalChromosomeContracts(
        final String rawOutdir,
        final Map references,
        final Map modalities
    ) {
        final File outputDirectory = new File(
            new File(rawOutdir).canonicalFile,
            'pipeline_info/derived_contract'
        )
        try {
            final Map contracts = [:]
            if( modalities.rna as boolean ) {
                contracts.rna = writeChromosomeContract(
                    'rna',
                    new File(references.rna_chrom_sizes.toString()),
                    'chrom-sizes',
                    outputDirectory
                )
            }
            if( modalities.dna as boolean ) {
                final String dnaChromSizes = references.dna_chrom_sizes?.toString()?.trim()
                contracts.dna = writeChromosomeContract(
                    'dna',
                    new File("${references.dna_bwa_reference}.ann"),
                    'bwa-ann',
                    outputDirectory,
                    dnaChromSizes ? new File(dnaChromSizes) : null
                )
            }
            if( contracts.isEmpty() ) {
                throw new IllegalArgumentException(
                    'At least one RNA or DNA reference dictionary is required'
                )
            }
            return contracts
        }
        catch( Exception error ) {
            throw new IllegalArgumentException(
                'Canonical chromosome resolution failed for the configured reference index: ' +
                error.message,
                error
            )
        }
    }

    static List<Map> readChromSizes(final File path) {
        final List<Map> entries = []
        final Set<String> seen = [] as Set<String>
        path.eachLine('UTF-8') { rawLine, lineNumber ->
            final String line = rawLine.trim()
            if( line && !line.startsWith('#') ) {
                final List<String> fields = line.split(/\s+/) as List<String>
                if( fields.size() < 2 ) {
                    throw new IllegalArgumentException(
                        "Malformed chromosome-size line ${lineNumber} in ${path}: '${rawLine}'"
                    )
                }
                long length
                try {
                    length = Long.parseLong(fields[1])
                }
                catch( NumberFormatException error ) {
                    throw new IllegalArgumentException(
                        "Invalid chromosome length on line ${lineNumber} in ${path}: '${fields[1]}'",
                        error
                    )
                }
                final String name = fields[0]
                if( length < 1 ) {
                    throw new IllegalArgumentException(
                        "Chromosome length must be positive on line ${lineNumber} in ${path}: ${length}"
                    )
                }
                if( !seen.add(name) ) {
                    throw new IllegalArgumentException("Duplicate chromosome '${name}' in ${path}")
                }
                entries << [name: name, length: length]
            }
        }
        if( entries.isEmpty() ) {
            throw new IllegalArgumentException("Reference chromosome dictionary is empty: ${path}")
        }
        return entries
    }

    static List<Map> readBwaAnn(final File path) {
        final List<String> lines = path.readLines('UTF-8')
            .collect { it.trim() }
            .findAll { it }
        if( lines.isEmpty() ) {
            throw new IllegalArgumentException("BWA annotation file is empty: ${path}")
        }

        final List<String> header = lines[0].split(/\s+/) as List<String>
        final int sequenceCount
        try {
            sequenceCount = Integer.parseInt(header[1])
        }
        catch( Exception error ) {
            throw new IllegalArgumentException(
                "Cannot read the sequence count from BWA annotation header in ${path}: '${lines[0]}'",
                error
            )
        }
        final int expectedLines = 1 + (sequenceCount * 2)
        if( sequenceCount < 1 || lines.size() < expectedLines ) {
            throw new IllegalArgumentException(
                "Malformed BWA annotation file ${path}: expected ${sequenceCount} sequence entries"
            )
        }

        final List<Map> entries = []
        final Set<String> seen = [] as Set<String>
        (0..<sequenceCount).each { index ->
            final List<String> descriptor = lines[1 + (index * 2)].split(/\s+/) as List<String>
            final List<String> coordinates = lines[2 + (index * 2)].split(/\s+/) as List<String>
            String name
            long length
            try {
                name = descriptor[1]
                length = Long.parseLong(coordinates[1])
            }
            catch( Exception error ) {
                throw new IllegalArgumentException(
                    "Malformed BWA sequence entry ${index + 1} in ${path}",
                    error
                )
            }
            if( length < 1 ) {
                throw new IllegalArgumentException(
                    "BWA chromosome length must be positive for '${name}' in ${path}: ${length}"
                )
            }
            if( !seen.add(name) ) {
                throw new IllegalArgumentException("Duplicate chromosome '${name}' in ${path}")
            }
            entries << [name: name, length: length]
        }
        return entries
    }

    static Map resolveCanonicalEntries(final List<Map> entries, final String sourceLabel) {
        final Set<String> names = entries.collect { it.name as String } as Set<String>
        final Set<String> ucscHits = names.intersect(UCSC_NAMES) as Set<String>
        final Set<String> ensemblHits = names.intersect(ENSEMBL_NAMES) as Set<String>

        if( ucscHits && ensemblHits ) {
            throw new IllegalArgumentException(
                'Cannot determine a coherent human chromosome naming convention for ' +
                "${sourceLabel}: both UCSC-style (${ucscHits.sort().join(', ')}) and " +
                "Ensembl-style (${ensemblHits.sort().join(', ')}) canonical names are present"
            )
        }
        if( !ucscHits && !ensemblHits ) {
            throw new IllegalArgumentException(
                "Cannot determine a human chromosome naming convention for ${sourceLabel}: " +
                'no exact canonical autosome, X, Y, or mitochondrial names were found'
            )
        }

        String style
        Set<String> allowedNames
        final List<String> missingAnchors = []
        if( ucscHits ) {
            style = 'ucsc'
            allowedNames = UCSC_NAMES
            if( !AUTOSOMES.any { chromosome -> names.contains('chr' + chromosome) } ) {
                missingAnchors << 'an autosome (chr1-chr22)'
            }
            if( !names.contains('chrX') ) missingAnchors << 'chrX'
            if( !names.contains('chrY') ) missingAnchors << 'chrY'
            if( !names.contains('chrM') ) missingAnchors << 'chrM'
        }
        else {
            style = 'ensembl'
            final List<String> mitochondrial = ['MT', 'M'].findAll { names.contains(it) }
            if( mitochondrial.size() > 1 ) {
                throw new IllegalArgumentException(
                    "Cannot choose a single Ensembl mitochondrial chromosome for ${sourceLabel}: " +
                    'both MT and M are present'
                )
            }
            allowedNames = (AUTOSOMES + ['X', 'Y'] + mitochondrial) as Set<String>
            if( !AUTOSOMES.any { chromosome -> names.contains(chromosome) } ) {
                missingAnchors << 'an autosome (1-22)'
            }
            if( !names.contains('X') ) missingAnchors << 'X'
            if( !names.contains('Y') ) missingAnchors << 'Y'
            if( mitochondrial.isEmpty() ) missingAnchors << 'MT or M'
        }

        if( missingAnchors ) {
            throw new IllegalArgumentException(
                "Cannot safely resolve canonical human chromosomes for ${sourceLabel}; missing: " +
                missingAnchors.join(', ')
            )
        }

        return [
            style  : style,
            entries: entries.findAll { allowedNames.contains(it.name as String) },
        ]
    }

    static void verifyMatchingContracts(
        final Map primary,
        final Map secondary,
        final String primaryLabel,
        final String secondaryLabel
    ) {
        if( primary.style != secondary.style ) {
            throw new IllegalArgumentException(
                "Reference naming mismatch: ${primaryLabel} is ${primary.style}-style but " +
                "${secondaryLabel} is ${secondary.style}-style"
            )
        }
        final Map primarySizes = primary.entries.collectEntries {
            [(it.name as String): it.length]
        }
        final Map secondarySizes = secondary.entries.collectEntries {
            [(it.name as String): it.length]
        }
        if( primarySizes != secondarySizes ) {
            throw new IllegalArgumentException(
                "Canonical chromosome names or lengths disagree between ${primaryLabel} and " +
                secondaryLabel
            )
        }
    }

    private static Map writeChromosomeContract(
        final String modality,
        final File source,
        final String sourceFormat,
        final File outputDirectory,
        final File verificationSizes = null
    ) {
        final List<Map> entries = sourceFormat == 'bwa-ann'
            ? readBwaAnn(source)
            : readChromSizes(source)
        final Map resolved = resolveCanonicalEntries(entries, source.toString())
        if( verificationSizes ) {
            final Map verification = resolveCanonicalEntries(
                readChromSizes(verificationSizes),
                verificationSizes.toString()
            )
            verifyMatchingContracts(
                resolved,
                verification,
                source.toString(),
                verificationSizes.toString()
            )
        }

        if( !outputDirectory.exists() && !outputDirectory.mkdirs() && !outputDirectory.exists() ) {
            throw new IllegalStateException(
                "Unable to create canonical chromosome output directory: ${outputDirectory}"
            )
        }
        final File allowlist = new File(outputDirectory, "${modality}_canonical_chromosomes.txt")
        final File chromSizes = new File(outputDirectory, "${modality}_canonical_chromosomes.chrom.sizes")
        allowlist.setText(
            resolved.entries.collect { "${it.name}\n" }.join(''),
            'UTF-8'
        )
        chromSizes.setText(
            resolved.entries.collect { "${it.name}\t${it.length}\n" }.join(''),
            'UTF-8'
        )
        return [
            style      : resolved.style,
            source     : source.canonicalPath,
            allowlist  : allowlist.canonicalPath,
            chrom_sizes: chromSizes.canonicalPath,
            contigs    : resolved.entries.collect { it.name },
        ]
    }

    static void writeRuntimeContract(final String rawOutdir) {
        final File pipelineInfoDir = new File((rawOutdir ?: 'results').toString(), 'pipeline_info')
        if( !pipelineInfoDir.exists() ) {
            pipelineInfoDir.mkdirs()
        }

        final File reportFile = new File(pipelineInfoDir, 'runtime_contract.tsv')
        reportFile.setText("tool\tconfigured_path\texists\tcurrently_used\n", 'UTF-8')
    }
}
