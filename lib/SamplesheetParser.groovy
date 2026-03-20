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
        if( !(parsed instanceof Map) || !parsed.samples ) {
            throw new IllegalArgumentException(
                "Samplesheet must contain a non-empty top-level 'samples:' mapping or list: ${samplesheetPath}"
            )
        }

        final File baseDir = sheetFile.parentFile ?: new File('.')
        final String libraryName = requireString(parsed.library_name, 'library_name')

        if( parsed.samples instanceof Map ) {
            return parseUnified(parsed, baseDir, libraryName, options)
        }

        if( parsed.samples instanceof List ) {
            return parseLegacy(parsed, baseDir, libraryName, options)
        }

        throw new IllegalArgumentException(
            "Top-level 'samples:' must be either a mapping of sample ids or a non-empty list"
        )
    }

    private static List<Map> parseUnified(
        final Map parsed,
        final File baseDir,
        final String libraryName,
        final Map options
    ) {
        final Map defaults = normalizedDefaults(options)
        final String ligationWhitelist = resolveExistingPath(
            baseDir,
            requireString(defaults.shared.ligation_whitelist, 'ligation_barcode_whitelist')
        )
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

    private static List<Map> parseLegacy(
        final Map parsed,
        final File baseDir,
        final String libraryName,
        final Map options
    ) {
        final List samples = (List) parsed.samples
        if( samples.isEmpty() ) {
            throw new IllegalArgumentException("Samplesheet must contain a non-empty top-level 'samples:' list")
        }

        final String sbGroupMap = resolveExistingPath(
            baseDir,
            requireString(parsed.sb_group_map, 'sb_group_map')
        )
        final boolean hasDnaSamples = samples.any { sample ->
            sample instanceof Map && sample.modality?.toString()?.trim()?.equalsIgnoreCase('dna')
        }
        final String dnaMoMap = hasDnaSamples
            ? resolveExistingPath(baseDir, requireString(parsed.dna_mo_map, 'dna_mo_map'))
            : null
        final Map defaults = normalizedDefaults(options)
        final String defaultLigationWhitelist = resolveExistingPath(
            baseDir,
            requireString(defaults.shared.ligation_whitelist, 'ligation_barcode_whitelist')
        )

        final List<Map> parsedRows = []

        for( int idx = 0; idx < samples.size(); idx++ ) {
            final def row = samples[idx]
            if( !(row instanceof Map) ) {
                throw new IllegalArgumentException("samples[${idx}] must be a mapping")
            }

            final String sampleId = requireString(row.id, "samples[${idx}].id")
            final String modality = requireString(row.modality, "samples[${idx}].modality").toLowerCase()

            final Map reads = asMap(row.reads, "samples[${idx}].reads")
            final Map barcodes = asMap(row.barcodes, "samples[${idx}].barcodes")
            final Map sampleBarcode = asMap(barcodes.sample, "samples[${idx}].barcodes.sample")
            final Map cellBarcode = barcodes.cell ? asMap(barcodes.cell, "samples[${idx}].barcodes.cell") : [:]

            if( modality == 'rna' ) {
                final Map umiBarcode = asMap(barcodes.umi, "samples[${idx}].barcodes.umi")

                parsedRows << [
                    id                         : sampleId,
                    modality                   : modality,
                    i1                         : resolveExistingPath(baseDir, requireString(reads.i1, "samples[${idx}].reads.i1")),
                    r1                         : resolveExistingPath(baseDir, requireString(reads.r1, "samples[${idx}].reads.r1")),
                    r2                         : resolveExistingPath(baseDir, requireString(reads.r2, "samples[${idx}].reads.r2")),
                    sample_bc_len              : requireInt(sampleBarcode.bc_len, "samples[${idx}].barcodes.sample.bc_len"),
                    sample_bc_start            : requireInt(sampleBarcode.bc_start, "samples[${idx}].barcodes.sample.bc_start"),
                    sample_hd                  : requireInt(sampleBarcode.hd, "samples[${idx}].barcodes.sample.hd"),
                    sample_tag                 : optionalString(sampleBarcode.tag, 'SB'),
                    sample_first_pass          : optionalString(sampleBarcode.first_pass, 'first_pass'),
                    sample_reverse_complement  : toDirection(
                        sampleBarcode.containsKey('reverse_complement') ? sampleBarcode.reverse_complement : true,
                        "samples[${idx}].barcodes.sample.reverse_complement"
                    ),
                    umi_bc_len                 : requireInt(umiBarcode.bc_len, "samples[${idx}].barcodes.umi.bc_len"),
                    umi_bc_start               : requireInt(umiBarcode.bc_start, "samples[${idx}].barcodes.umi.bc_start"),
                    umi_tag                    : optionalString(umiBarcode.tag, 'UM'),
                    cell_whitelist             : resolveLegacyCellWhitelist(baseDir, cellBarcode, defaultLigationWhitelist, idx),
                    cell_bc_len                : requireInt(cellBarcode.bc_len ?: defaults.rna.cell.bc_len, "samples[${idx}].barcodes.cell.bc_len"),
                    cell_hd                    : requireInt(cellBarcode.hd ?: defaults.rna.cell.hd, "samples[${idx}].barcodes.cell.hd"),
                    cell_tag                   : optionalString(cellBarcode.tag, defaults.rna.cell.tag.toString()),
                    library_name               : libraryName,
                    sb_group_map               : sbGroupMap
                ]
                continue
            }

            if( modality == 'dna' ) {
                final Map modalityBarcode = asMap(barcodes.modality, "samples[${idx}].barcodes.modality")

                parsedRows << [
                    id                            : sampleId,
                    modality                      : modality,
                    i1                            : resolveExistingPath(baseDir, requireString(reads.i1, "samples[${idx}].reads.i1")),
                    i2                            : resolveExistingPath(baseDir, requireString(reads.i2, "samples[${idx}].reads.i2")),
                    r1                            : resolveExistingPath(baseDir, requireString(reads.r1, "samples[${idx}].reads.r1")),
                    r2                            : resolveExistingPath(baseDir, requireString(reads.r2, "samples[${idx}].reads.r2")),
                    sample_bc_len                 : requireInt(sampleBarcode.bc_len, "samples[${idx}].barcodes.sample.bc_len"),
                    sample_bc_start               : requireInt(sampleBarcode.bc_start, "samples[${idx}].barcodes.sample.bc_start"),
                    sample_hd                     : requireInt(sampleBarcode.hd, "samples[${idx}].barcodes.sample.hd"),
                    sample_tag                    : optionalString(sampleBarcode.tag, 'SB'),
                    sample_first_pass             : optionalString(sampleBarcode.first_pass, 'first_pass'),
                    sample_reverse_complement     : toDirection(
                        sampleBarcode.containsKey('reverse_complement') ? sampleBarcode.reverse_complement : true,
                        "samples[${idx}].barcodes.sample.reverse_complement"
                    ),
                    modality_whitelist            : resolveExistingPath(
                        baseDir,
                        requireString(modalityBarcode.whitelist, "samples[${idx}].barcodes.modality.whitelist")
                    ),
                    modality_bc_len               : requireInt(modalityBarcode.bc_len, "samples[${idx}].barcodes.modality.bc_len"),
                    modality_bc_start             : requireInt(modalityBarcode.bc_start, "samples[${idx}].barcodes.modality.bc_start"),
                    modality_hd                   : requireInt(modalityBarcode.hd, "samples[${idx}].barcodes.modality.hd"),
                    modality_tag                  : optionalString(modalityBarcode.tag, 'MO'),
                    modality_first_pass           : optionalString(modalityBarcode.first_pass, 'not_first_pass'),
                    modality_reverse_complement   : toDirection(
                        modalityBarcode.containsKey('reverse_complement') ? modalityBarcode.reverse_complement : true,
                        "samples[${idx}].barcodes.modality.reverse_complement"
                    ),
                    cell_whitelist                : resolveLegacyCellWhitelist(baseDir, cellBarcode, defaultLigationWhitelist, idx),
                    cell_bc_len                   : requireInt(cellBarcode.bc_len ?: defaults.dna.cell.bc_len, "samples[${idx}].barcodes.cell.bc_len"),
                    cell_hd                       : requireInt(cellBarcode.hd ?: defaults.dna.cell.hd, "samples[${idx}].barcodes.cell.hd"),
                    cell_tag                      : optionalString(cellBarcode.tag, defaults.dna.cell.tag.toString()),
                    library_name                  : libraryName,
                    sb_group_map                  : sbGroupMap,
                    mo_map                        : dnaMoMap
                ]
                continue
            }

            throw new IllegalArgumentException(
                "Unsupported modality for sample '${sampleId}': '${modality}'. Supported values: rna, dna"
            )
        }

        return parsedRows
    }

    private static String resolveLegacyCellWhitelist(
        final File baseDir,
        final Map cellBarcode,
        final String defaultLigationWhitelist,
        final int idx
    ) {
        if( cellBarcode?.whitelist ) {
            return resolveExistingPath(
                baseDir,
                requireString(cellBarcode.whitelist, "samples[${idx}].barcodes.cell.whitelist")
            )
        }
        return defaultLigationWhitelist
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
            shared: [
                ligation_whitelist: options.ligation_barcode_whitelist
            ],
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

    private static int requireInt(final Object value, final String fieldName) {
        if( value == null ) {
            throw new IllegalArgumentException("Missing required field: ${fieldName}")
        }
        try {
            return value as int
        }
        catch( Exception ignored ) {
            throw new IllegalArgumentException("${fieldName} must be an integer")
        }
    }

    private static String optionalString(final Object value, final String defaultValue) {
        final String out = value?.toString()?.trim()
        return out ?: defaultValue
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
}
