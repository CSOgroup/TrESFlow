import groovy.yaml.YamlSlurper

class SamplesheetParser {

    static List<Map> parse(final String samplesheetPath, final Map options = [:]) {
        if( !samplesheetPath ) {
            throw new IllegalArgumentException("Missing required parameter: --samplesheet")
        }

        final File sheetFile = new File(samplesheetPath)
        if( !sheetFile.exists() ) {
            throw new IllegalArgumentException("Samplesheet not found: ${samplesheetPath}")
        }

        final def parsed = new YamlSlurper().parse(sheetFile)
        if( !(parsed instanceof Map) ) {
            throw new IllegalArgumentException("Samplesheet must be a top-level YAML mapping: ${samplesheetPath}")
        }
        if( !(parsed.resources instanceof Map) || ((Map) parsed.resources).isEmpty() ) {
            throw new IllegalArgumentException(
                "Samplesheet must contain a non-empty top-level 'resources:' mapping: ${samplesheetPath}"
            )
        }
        if( !(parsed.samples instanceof Map) || ((Map) parsed.samples).isEmpty() ) {
            throw new IllegalArgumentException(
                "Samplesheet must contain a non-empty top-level 'samples:' mapping: ${samplesheetPath}"
            )
        }

        final File baseDir = sheetFile.parentFile ?: new File('.')
        final String libraryName = requireString(parsed.library_name, 'library_name')
        final Map sharedResources = resolveSharedResources(parsed, baseDir, options)

        return parseUnified(parsed, baseDir, libraryName, options, sharedResources)
    }

    private static List<Map> parseUnified(
        final Map parsed,
        final File baseDir,
        final String libraryName,
        final Map options,
        final Map sharedResources
    ) {
        final Map defaults = normalizedDefaults(options)
        final String ligationWhitelist = sharedResources.ligation_barcode_whitelist
        final File derivedDir = prepareDerivedDir(options)

        final List<Map> samples = []

        ((Map) parsed.samples).each { rawSampleId, rawSampleConfig ->
            final String sampleId = requireString(rawSampleId, 'samples.<sample_id>')
            final Map sampleConfig = asMap(rawSampleConfig, "samples.${sampleId}")
            final Map groupsConfig = asMap(sampleConfig.groups, "samples.${sampleId}.groups")
            if( groupsConfig.isEmpty() ) {
                throw new IllegalArgumentException("samples.${sampleId}.groups must not be empty")
            }

            final LinkedHashMap<String, List<String>> normalizedGroups = parseGroups(groupsConfig, sampleId)

            final Map rnaConfig = sampleConfig.rna ? asMap(sampleConfig.rna, "samples.${sampleId}.rna") : null
            final Map dnaConfig = sampleConfig.dna ? asMap(sampleConfig.dna, "samples.${sampleId}.dna") : null
            if( !rnaConfig && !dnaConfig ) {
                throw new IllegalArgumentException(
                    "samples.${sampleId} must define at least one modality block: rna or dna"
                )
            }

            if( rnaConfig ) {
                final Map reads = asMap(rnaConfig.reads, "samples.${sampleId}.rna.reads")
                samples << [
                    id                        : sampleId,
                    modality                  : 'rna',
                    i1                        : resolveExistingPath(baseDir, requireString(reads.i1, "samples.${sampleId}.rna.reads.i1")),
                    r1                        : resolveExistingPath(baseDir, requireString(reads.r1, "samples.${sampleId}.rna.reads.r1")),
                    r2                        : resolveExistingPath(baseDir, requireString(reads.r2, "samples.${sampleId}.rna.reads.r2")),
                    sample_bc_len             : defaults.rna.sample.bc_len as int,
                    sample_bc_start           : defaults.rna.sample.bc_start as int,
                    sample_hd                 : defaults.rna.sample.hd as int,
                    sample_tag                : defaults.rna.sample.tag.toString(),
                    sample_first_pass         : defaults.rna.sample.first_pass.toString(),
                    sample_reverse_complement : toDirection(defaults.rna.sample.reverse_complement, 'barcode_defaults.rna.sample.reverse_complement'),
                    umi_bc_len                : defaults.rna.umi.bc_len as int,
                    umi_bc_start              : defaults.rna.umi.bc_start as int,
                    umi_tag                   : defaults.rna.umi.tag.toString(),
                    cell_whitelist            : ligationWhitelist,
                    cell_bc_len               : defaults.rna.cell.bc_len as int,
                    cell_hd                   : defaults.rna.cell.hd as int,
                    cell_tag                  : defaults.rna.cell.tag.toString(),
                    rna_ref_base_dir          : sharedResources.rna_ref_base_dir,
                    rna_align_species         : sharedResources.rna_align_species,
                    library_name              : libraryName,
                    group_definitions         : normalizedGroups
                ]
            }

            if( dnaConfig ) {
                final Map reads = asMap(dnaConfig.reads, "samples.${sampleId}.dna.reads")
                final LinkedHashMap<String, String> markBarcodes = parseMarkBarcodes(
                    dnaConfig.mark_barcodes,
                    sampleId
                )

                samples << [
                    id                          : sampleId,
                    modality                    : 'dna',
                    i1                          : resolveExistingPath(baseDir, requireString(reads.i1, "samples.${sampleId}.dna.reads.i1")),
                    i2                          : resolveExistingPath(baseDir, requireString(reads.i2, "samples.${sampleId}.dna.reads.i2")),
                    r1                          : resolveExistingPath(baseDir, requireString(reads.r1, "samples.${sampleId}.dna.reads.r1")),
                    r2                          : resolveExistingPath(baseDir, requireString(reads.r2, "samples.${sampleId}.dna.reads.r2")),
                    sample_bc_len               : defaults.dna.sample.bc_len as int,
                    sample_bc_start             : defaults.dna.sample.bc_start as int,
                    sample_hd                   : defaults.dna.sample.hd as int,
                    sample_tag                  : defaults.dna.sample.tag.toString(),
                    sample_first_pass           : defaults.dna.sample.first_pass.toString(),
                    sample_reverse_complement   : toDirection(defaults.dna.sample.reverse_complement, 'barcode_defaults.dna.sample.reverse_complement'),
                    modality_bc_len             : defaults.dna.modality.bc_len as int,
                    modality_bc_start           : defaults.dna.modality.bc_start as int,
                    modality_hd                 : defaults.dna.modality.hd as int,
                    modality_tag                : defaults.dna.modality.tag.toString(),
                    modality_first_pass         : defaults.dna.modality.first_pass.toString(),
                    modality_reverse_complement : toDirection(defaults.dna.modality.reverse_complement, 'barcode_defaults.dna.modality.reverse_complement'),
                    cell_whitelist              : ligationWhitelist,
                    cell_bc_len                 : defaults.dna.cell.bc_len as int,
                    cell_hd                     : defaults.dna.cell.hd as int,
                    cell_tag                    : defaults.dna.cell.tag.toString(),
                    dna_bwa_reference           : sharedResources.dna_bwa_reference,
                    dna_blacklist_bed           : sharedResources.dna_blacklist_bed,
                    dna_effective_genome_size   : sharedResources.dna_effective_genome_size,
                    library_name                : libraryName,
                    group_definitions           : normalizedGroups,
                    mark_barcodes               : markBarcodes
                ]
            }
        }

        final File sbGroupMapFile = writeSbGroupMap(derivedDir, samples)
        final List<Map> dnaRows = samples.findAll { it.modality == 'dna' }
        final File dnaMoMapFile = dnaRows ? writeDnaMoMap(derivedDir, dnaRows) : null
        final Map<String, String> dnaWhitelistPaths = dnaRows ? writeDnaModalityWhitelists(derivedDir, dnaRows) : [:]

        samples.each { row ->
            row.sb_group_map = sbGroupMapFile.canonicalPath
            row.remove('group_definitions')
            if( row.modality == 'dna' ) {
                row.mo_map = dnaMoMapFile.canonicalPath
                row.modality_whitelist = dnaWhitelistPaths[row.id]
                row.remove('mark_barcodes')
            }
        }

        return samples
    }

    private static LinkedHashMap<String, List<String>> parseGroups(final Map groupsConfig, final String sampleId) {
        final LinkedHashMap<String, List<String>> normalized = new LinkedHashMap<>()
        final Set<String> seenBarcodes = new LinkedHashSet<>()

        groupsConfig.each { rawGroupName, rawGroupConfig ->
            final String groupName = requireString(rawGroupName, "samples.${sampleId}.groups.<group>")
            final Map groupConfig = asMap(rawGroupConfig, "samples.${sampleId}.groups.${groupName}")
            final List<String> barcodes = requireBarcodeList(
                groupConfig.sb_barcodes,
                "samples.${sampleId}.groups.${groupName}.sb_barcodes"
            )

            barcodes.each { barcode ->
                if( !seenBarcodes.add(barcode) ) {
                    throw new IllegalArgumentException(
                        "Duplicate sample barcode '${barcode}' across groups for sample '${sampleId}'"
                    )
                }
            }

            normalized[groupName] = barcodes
        }

        return normalized
    }

    private static LinkedHashMap<String, String> parseMarkBarcodes(final Object value, final String sampleId) {
        final Map marksConfig = asMap(value, "samples.${sampleId}.dna.mark_barcodes")
        if( marksConfig.isEmpty() ) {
            throw new IllegalArgumentException("samples.${sampleId}.dna.mark_barcodes must not be empty")
        }

        final LinkedHashMap<String, String> normalized = new LinkedHashMap<>()
        final Map<String, String> barcodeToMark = [:]

        marksConfig.each { rawMarkName, rawBarcode ->
            final String markName = requireString(rawMarkName, "samples.${sampleId}.dna.mark_barcodes.<mark>")
            final String barcode = requireString(rawBarcode, "samples.${sampleId}.dna.mark_barcodes.${markName}")
            if( barcodeToMark.containsKey(barcode) && barcodeToMark[barcode] != markName ) {
                throw new IllegalArgumentException(
                    "Duplicate DNA modality barcode '${barcode}' for sample '${sampleId}': " +
                    "${barcodeToMark[barcode]} vs ${markName}"
                )
            }
            barcodeToMark[barcode] = markName
            normalized[markName] = barcode
        }

        return normalized
    }

    private static List<String> requireBarcodeList(final Object value, final String fieldName) {
        if( !(value instanceof List) || value.isEmpty() ) {
            throw new IllegalArgumentException("${fieldName} must be a non-empty list")
        }

        final List<String> normalized = value.collect { entry ->
            requireString(entry, fieldName)
        }

        if( normalized.toSet().size() != normalized.size() ) {
            throw new IllegalArgumentException("${fieldName} contains duplicate barcodes")
        }

        return normalized
    }

    private static File prepareDerivedDir(final Map options) {
        final String outdir = options.outdir?.toString()?.trim()
        final File root = outdir
            ? new File(new File(outdir), 'pipeline_info/derived_contract')
            : File.createTempDir('tresflow_samplesheet_', '')

        if( root.exists() ) {
            root.eachFile { file -> deleteRecursively(file) }
        }
        root.mkdirs()
        return root
    }

    private static void deleteRecursively(final File file) {
        if( file.isDirectory() ) {
            file.listFiles()?.each { child -> deleteRecursively(child) }
        }
        file.delete()
    }

    private static File writeSbGroupMap(final File derivedDir, final List<Map> samples) {
        final File file = new File(derivedDir, 'sb_group_map.tsv')
        final List<String> lines = ['sample\tsb_group\tsb_bc']

        samples.collect { it.id }.unique().each { sampleId ->
            final Map row = samples.find { it.id == sampleId }
            row.group_definitions.each { groupName, barcodes ->
                barcodes.each { barcode ->
                    lines << "${sampleId}\t${groupName}\t${barcode}"
                }
            }
        }

        file.text = lines.join('\n') + '\n'
        return file
    }

    private static File writeDnaMoMap(final File derivedDir, final List<Map> dnaRows) {
        final File file = new File(derivedDir, 'dna_mo_map.tsv')
        final List<String> lines = ['sample\tsb_group\tmark\tmo_bc']

        dnaRows.each { row ->
            row.group_definitions.each { groupName, barcodes ->
                row.mark_barcodes.each { markName, barcode ->
                    lines << "${row.id}\t${groupName}\t${markName}\t${barcode}"
                }
            }
        }

        file.text = lines.join('\n') + '\n'
        return file
    }

    private static Map<String, String> writeDnaModalityWhitelists(final File derivedDir, final List<Map> dnaRows) {
        final File whitelistDir = new File(derivedDir, 'dna_modality_whitelists')
        whitelistDir.mkdirs()

        final Map<String, String> out = [:]
        dnaRows.each { row ->
            final File file = new File(whitelistDir, "${row.id}.txt")
            file.text = row.mark_barcodes.values().join('\n') + '\n'
            out[row.id] = file.canonicalPath
        }
        return out
    }

    private static Map normalizedDefaults(final Map options) {
        final Map barcodeDefaults = asMap(options.barcode_defaults ?: [:], 'barcode_defaults')
        return [
            rna: [
                sample: asMap(barcodeDefaults.rna?.sample, 'barcode_defaults.rna.sample'),
                umi   : asMap(barcodeDefaults.rna?.umi, 'barcode_defaults.rna.umi'),
                cell  : asMap(barcodeDefaults.rna?.cell, 'barcode_defaults.rna.cell')
            ],
            dna: [
                sample  : asMap(barcodeDefaults.dna?.sample, 'barcode_defaults.dna.sample'),
                modality: asMap(barcodeDefaults.dna?.modality, 'barcode_defaults.dna.modality'),
                cell    : asMap(barcodeDefaults.dna?.cell, 'barcode_defaults.dna.cell')
            ]
        ]
    }

    private static Map resolveSharedResources(final Map parsed, final File baseDir, final Map options) {
        final Map resources = asMap(parsed.resources, 'resources')

        return [
            ligation_barcode_whitelist: resolveExistingPath(
                baseDir,
                requireString(
                    firstDefined(resources.ligation_barcode_whitelist, options.ligation_barcode_whitelist),
                    'resources.ligation_barcode_whitelist or --ligation_barcode_whitelist'
                )
            ),
            rna_ref_base_dir: resolveOptionalPath(
                baseDir,
                firstDefined(resources.rna_ref_base_dir, options.rna_ref_base_dir)
            ),
            rna_align_species: firstDefined(
                resources.rna_align_species,
                options.rna_align_species
            )?.toString()?.trim()?.toLowerCase(),
            dna_bwa_reference: resolveOptionalPath(
                baseDir,
                firstDefined(resources.dna_bwa_reference, options.dna_bwa_reference)
            ),
            dna_blacklist_bed: resolveOptionalPath(
                baseDir,
                firstDefined(resources.dna_blacklist_bed, options.dna_blacklist_bed)
            ),
            dna_effective_genome_size: firstDefined(
                resources.dna_effective_genome_size,
                options.dna_effective_genome_size
            )?.toString()?.trim()
        ]
    }

    private static Object firstDefined(final Object... values) {
        for( final Object value : values ) {
            if( value == null ) {
                continue
            }
            if( value instanceof CharSequence && !value.toString().trim() ) {
                continue
            }
            return value
        }
        return null
    }

    private static Map asMap(final Object value, final String fieldName) {
        if( value instanceof Map ) {
            return (Map) value
        }
        throw new IllegalArgumentException("${fieldName} must be a mapping")
    }

    private static String requireString(final Object value, final String fieldName) {
        final String out = value?.toString()?.trim()
        if( !out ) {
            throw new IllegalArgumentException("Missing required field: ${fieldName}")
        }
        return out
    }

    private static String toDirection(final Object value, final String fieldName) {
        if( value instanceof Boolean ) {
            return value ? 'rev' : 'fw'
        }

        final String normalized = value?.toString()?.trim()?.toLowerCase()
        if( normalized in ['true', 'rev', 'reverse', 'reverse_complement'] ) {
            return 'rev'
        }
        if( normalized in ['false', 'fw', 'forward'] ) {
            return 'fw'
        }

        throw new IllegalArgumentException(
            "${fieldName} must be a boolean or one of: rev, fw, reverse, forward"
        )
    }

    private static String resolveExistingPath(final File baseDir, final String rawPath) {
        final File resolved = new File(rawPath).isAbsolute() ? new File(rawPath) : new File(baseDir, rawPath)
        if( !resolved.exists() ) {
            throw new IllegalArgumentException("Referenced file not found: ${resolved}")
        }
        return resolved.canonicalPath
    }

    private static String resolveOptionalPath(final File baseDir, final Object rawPath) {
        final String normalized = rawPath?.toString()?.trim()
        if( !normalized ) {
            return null
        }

        final File resolved = new File(normalized).isAbsolute() ? new File(normalized) : new File(baseDir, normalized)
        return resolved.canonicalPath
    }
}
