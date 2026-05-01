import groovy.yaml.YamlSlurper

class SamplesheetParser {

    private static final String MODALITY_RNA = 'rna'
    private static final String MODALITY_DNA = 'dna'
    private static final String TAGMENTATION_SINGLE = 'single'
    private static final String TAGMENTATION_DUAL = 'dual'
    private static final List<String> DNA_TAGMENTATION_MODES = [TAGMENTATION_SINGLE, TAGMENTATION_DUAL]
    private static final int RNA_SB_BARCODE_LENGTH = 4
    private static final int DNA_SINGLE_SB_BARCODE_LENGTH = 4
    private static final int DNA_DUAL_SB_BARCODE_LENGTH = 3

    static List<Map> parse(final String samplesheetPath, final Map options = [:]) {
        return parseContract(samplesheetPath, options).samples as List<Map>
    }

    static Map parseContract(final String samplesheetPath, final Map options = [:]) {
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
        if( !(parsed.samples instanceof Map) || ((Map) parsed.samples).isEmpty() ) {
            throw new IllegalArgumentException(
                "Samplesheet must contain a non-empty top-level 'samples:' mapping: ${samplesheetPath}"
            )
        }

        final File baseDir = sheetFile.parentFile ?: new File('.')
        final String libraryName = requireString(parsed.library_name, 'library_name')
        final Map runtime = resolveRuntime(parsed, baseDir)
        final Map references = resolveReferences(parsed, baseDir)
        final List<Map> samples = parseUnified(parsed, baseDir, libraryName, options, references)
        samples.each { row ->
            row.runtime_env_prefix = runtime['env_prefix']
            row.runtime_tmpdir = runtime['tmpdir']
        }

        return [
            library_name: libraryName,
            runtime     : runtime,
            references  : references,
            modalities  : [
                rna: samples.any { row -> row.modality == MODALITY_RNA },
                dna: samples.any { row -> row.modality == MODALITY_DNA },
            ],
            samples     : samples,
        ]
    }

    private static List<Map> parseUnified(
        final Map parsed,
        final File baseDir,
        final String libraryName,
        final Map options,
        final Map references
    ) {
        // Internally the workflow still runs one modality row at a time. The
        // public contract stays hierarchical; the parser is where that view is
        // flattened into RNA and DNA work rows plus derived helper files.
        final Map defaults = normalizedDefaults(options)
        final String ligationWhitelist = references.ligation_barcode_whitelist
        final File derivedDir = prepareDerivedDir(options)

        final List<Map> samples = []

        ((Map) parsed.samples).each { rawSampleId, rawSampleConfig ->
            final String sampleId = requireString(rawSampleId, 'samples.<sample_id>')
            final Map sampleConfig = asMap(rawSampleConfig, "samples.${sampleId}")
            final Map groupsConfig = asMap(sampleConfig.groups, "samples.${sampleId}.groups")
            if( groupsConfig.isEmpty() ) {
                throw new IllegalArgumentException("samples.${sampleId}.groups must not be empty")
            }

            final Map rnaConfig = sampleConfig.rna ? asMap(sampleConfig.rna, "samples.${sampleId}.rna") : null
            final Map dnaConfig = sampleConfig.dna ? asMap(sampleConfig.dna, "samples.${sampleId}.dna") : null
            final boolean hasRna = rnaConfig != null
            final boolean hasDna = dnaConfig != null

            if( !hasRna && !hasDna ) {
                throw new IllegalArgumentException(
                    "samples.${sampleId} must define at least one modality block: rna or dna"
                )
            }
            final String dnaTagmentation = hasDna ? parseDnaTagmentation(dnaConfig, sampleId) : null
            final Map normalizedGroups = parseGroups(
                groupsConfig,
                sampleId,
                hasRna,
                hasDna,
                dnaTagmentation
            )

            if( hasRna ) {
                samples << buildRnaRow(
                    sampleId,
                    libraryName,
                    baseDir,
                    normalizedGroups.rna as LinkedHashMap<String, List<String>>,
                    normalizedGroups.rna_source_summary as String,
                    rnaConfig,
                    defaults,
                    ligationWhitelist,
                    references
                )
            }

            if( hasDna ) {
                samples << buildDnaRow(
                    sampleId,
                    libraryName,
                    baseDir,
                    normalizedGroups.dna as LinkedHashMap<String, List<String>>,
                    normalizedGroups.dna_source_summary as String,
                    dnaTagmentation,
                    dnaConfig,
                    defaults,
                    ligationWhitelist,
                    references
                )
            }
        }

        return attachDerivedArtifacts(derivedDir, samples)
    }

    private static Map buildRnaRow(
        final String sampleId,
        final String libraryName,
        final File baseDir,
        final LinkedHashMap<String, List<String>> normalizedGroups,
        final String sbSourceSummary,
        final Map rnaConfig,
        final Map defaults,
        final String ligationWhitelist,
        final Map references
    ) {
        final Map reads = resolveReads(baseDir, rnaConfig, sampleId, MODALITY_RNA, ['i1', 'r1', 'r2'])
        return [
            id                        : sampleId,
            modality                  : MODALITY_RNA,
            i1                        : reads.i1,
            r1                        : reads.r1,
            r2                        : reads.r2,
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
            reference_species         : references.species,
            rna_star_index_dir        : requireString(references.rna_ref_dir, 'references.rna_ref_dir'),
            rna_chrom_sizes           : requireString(references.rna_chrom_sizes, 'references.rna_ref_dir/chrNameLength.txt'),
            library_name              : libraryName,
            group_definitions         : normalizedGroups,
            rna_sb_barcode_source     : sbSourceSummary,
            rna_sb_barcode_len        : RNA_SB_BARCODE_LENGTH,
        ]
    }

    private static Map buildDnaRow(
        final String sampleId,
        final String libraryName,
        final File baseDir,
        final LinkedHashMap<String, List<String>> normalizedGroups,
        final String sbSourceSummary,
        final String tagmentation,
        final Map dnaConfig,
        final Map defaults,
        final String ligationWhitelist,
        final Map references
    ) {
        final Map reads = resolveReads(baseDir, dnaConfig, sampleId, MODALITY_DNA, ['i1', 'i2', 'r1', 'r2'])
        final LinkedHashMap<String, String> markBarcodes = parseMarkBarcodes(
            dnaConfig.mark_barcodes,
            sampleId
        )
        final Map dnaTagDefaults = dnaTagDefaultsForTagmentation(defaults, tagmentation)

        return [
            id                          : sampleId,
            modality                    : MODALITY_DNA,
            dna_tagmentation            : tagmentation,
            i1                          : reads.i1,
            i2                          : reads.i2,
            r1                          : reads.r1,
            r2                          : reads.r2,
            sample_bc_len               : dnaTagDefaults.sample.bc_len as int,
            sample_bc_start             : dnaTagDefaults.sample.bc_start as int,
            sample_hd                   : dnaTagDefaults.sample.hd as int,
            sample_tag                  : defaults.dna.sample.tag.toString(),
            sample_first_pass           : defaults.dna.sample.first_pass.toString(),
            sample_reverse_complement   : dnaTagDefaults.sample.reverse_complement.toString(),
            modality_bc_len             : dnaTagDefaults.modality.bc_len as int,
            modality_bc_start           : dnaTagDefaults.modality.bc_start as int,
            modality_hd                 : dnaTagDefaults.modality.hd as int,
            modality_tag                : defaults.dna.modality.tag.toString(),
            modality_first_pass         : defaults.dna.modality.first_pass.toString(),
            modality_reverse_complement : dnaTagDefaults.modality.reverse_complement.toString(),
            dna_sample_index_read       : dnaTagDefaults.sample.index_read.toString(),
            dna_modality_index_read     : dnaTagDefaults.modality.index_read.toString(),
            cell_whitelist              : ligationWhitelist,
            cell_bc_len                 : defaults.dna.cell.bc_len as int,
            cell_hd                     : defaults.dna.cell.hd as int,
            cell_tag                    : defaults.dna.cell.tag.toString(),
            reference_species           : references.species,
            dna_ref_dir                 : references.dna_ref_dir,
            dna_bwa_reference           : references.dna_bwa_reference ?: '',
            dna_blacklist_bed           : references.dna_blacklist_bed,
            dna_chrom_sizes             : references.dna_chrom_sizes,
            dna_effective_genome_size   : references.dna_effective_genome_size,
            library_name                : libraryName,
            group_definitions           : normalizedGroups,
            dna_sb_barcode_source       : sbSourceSummary,
            dna_sb_barcode_len          : dnaSbBarcodeLength(tagmentation),
            mark_barcodes               : markBarcodes,
        ]
    }

    private static List<Map> attachDerivedArtifacts(final File derivedDir, final List<Map> samples) {
        final List<Map> rnaRows = samples.findAll { it.modality == MODALITY_RNA }
        final List<Map> dnaRows = samples.findAll { it.modality == MODALITY_DNA }
        final File rnaSbGroupMapFile = rnaRows ? writeSbGroupMap(derivedDir, rnaRows, MODALITY_RNA) : null
        final File dnaSbGroupMapFile = dnaRows ? writeSbGroupMap(derivedDir, dnaRows, MODALITY_DNA) : null
        final File dnaMoMapFile = dnaRows ? writeDnaMoMap(derivedDir, dnaRows) : null
        final Map<String, String> dnaWhitelistPaths = dnaRows ? writeDnaModalityWhitelists(derivedDir, dnaRows) : [:]

        samples.each { row ->
            row.sb_group_map = row.modality == MODALITY_RNA
                ? rnaSbGroupMapFile.canonicalPath
                : dnaSbGroupMapFile.canonicalPath
            row.remove('group_definitions')
            if( row.modality == MODALITY_DNA ) {
                row.mo_map = dnaMoMapFile.canonicalPath
                row.modality_whitelist = dnaWhitelistPaths[row.id]
                row.remove('mark_barcodes')
            }
        }

        return samples
    }

    private static Map parseGroups(
        final Map groupsConfig,
        final String sampleId,
        final boolean hasRna,
        final boolean hasDna,
        final String dnaTagmentation
    ) {
        final LinkedHashMap<String, List<String>> rnaGroups = new LinkedHashMap<>()
        final LinkedHashMap<String, List<String>> dnaGroups = new LinkedHashMap<>()
        final LinkedHashMap<String, String> rnaSources = new LinkedHashMap<>()
        final LinkedHashMap<String, String> dnaSources = new LinkedHashMap<>()

        groupsConfig.each { rawGroupName, rawGroupConfig ->
            final String groupName = requireString(rawGroupName, "samples.${sampleId}.groups.<group>")
            final Map groupConfig = asMap(rawGroupConfig, "samples.${sampleId}.groups.${groupName}")

            if( hasRna ) {
                final Map rnaSelection = selectRnaSbBarcodes(groupConfig, sampleId, groupName)
                addGroupBarcodes(rnaGroups, rnaSources, groupName, rnaSelection, RNA_SB_BARCODE_LENGTH)
            }

            if( hasDna ) {
                final Map dnaSelection = selectDnaSbBarcodes(groupConfig, sampleId, groupName, dnaTagmentation)
                addGroupBarcodes(dnaGroups, dnaSources, groupName, dnaSelection, dnaSbBarcodeLength(dnaTagmentation))
            }
        }

        if( hasRna ) {
            validateNoGroupBarcodeCollisions(sampleId, 'RNA', rnaGroups)
        }
        if( hasDna ) {
            validateNoGroupBarcodeCollisions(sampleId, 'DNA', dnaGroups)
        }

        return [
            rna               : rnaGroups,
            dna               : dnaGroups,
            rna_source_summary: summarizeSources(rnaSources.values()),
            dna_source_summary: summarizeSources(dnaSources.values()),
        ]
    }

    private static void addGroupBarcodes(
        final LinkedHashMap<String, List<String>> groups,
        final LinkedHashMap<String, String> sources,
        final String groupName,
        final Map selection,
        final int expectedLength
    ) {
        final String fieldName = selection.fieldName as String
        groups[groupName] = requireBarcodeList(selection.value, fieldName, expectedLength)
        sources[groupName] = fieldName
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

    private static Map selectRnaSbBarcodes(final Map groupConfig, final String sampleId, final String groupName) {
        final String explicitField = "samples.${sampleId}.groups.${groupName}.rna_sb_barcodes"
        if( groupConfig.containsKey('rna_sb_barcodes') ) {
            return [value: groupConfig.rna_sb_barcodes, fieldName: explicitField]
        }

        return [
            value    : groupConfig.sb_barcodes,
            fieldName: "samples.${sampleId}.groups.${groupName}.sb_barcodes",
        ]
    }

    private static Map selectDnaSbBarcodes(
        final Map groupConfig,
        final String sampleId,
        final String groupName,
        final String dnaTagmentation
    ) {
        final String explicitField = "samples.${sampleId}.groups.${groupName}.dna_sb_barcodes"
        if( groupConfig.containsKey('dna_sb_barcodes') ) {
            return [value: groupConfig.dna_sb_barcodes, fieldName: explicitField]
        }

        if( dnaTagmentation == TAGMENTATION_DUAL ) {
            throw new IllegalArgumentException(
                "Missing required field: ${explicitField}. DNA dual tagmentation requires explicit 3 nt dna_sb_barcodes; " +
                "the pipeline will not derive them from RNA/sample sb_barcodes."
            )
        }

        return [
            value    : groupConfig.sb_barcodes,
            fieldName: "samples.${sampleId}.groups.${groupName}.sb_barcodes",
        ]
    }

    private static void validateNoGroupBarcodeCollisions(
        final String sampleId,
        final String modality,
        final LinkedHashMap<String, List<String>> groups
    ) {
        final Map<String, String> barcodeToGroup = [:]
        groups.each { groupName, barcodes ->
            barcodes.each { barcode ->
                if( barcodeToGroup.containsKey(barcode) && barcodeToGroup[barcode] != groupName ) {
                    throw new IllegalArgumentException(
                        "${modality} sample barcode collision for sample '${sampleId}': barcode '${barcode}' " +
                        "maps to both group '${barcodeToGroup[barcode]}' and group '${groupName}'"
                    )
                }
                barcodeToGroup[barcode] = groupName
            }
        }
    }

    private static String summarizeSources(final Collection<String> sources) {
        return sources.findAll { it }.collect { source ->
            source.tokenize('.').last()
        }.unique().sort().join(',')
    }

    private static List<String> requireBarcodeList(final Object value, final String fieldName, final int expectedLength) {
        if( !(value instanceof List) || value.isEmpty() ) {
            throw new IllegalArgumentException("${fieldName} must be a non-empty list")
        }

        final List<String> normalized = value.collect { entry ->
            requireString(entry, fieldName)
        }

        if( normalized.toSet().size() != normalized.size() ) {
            throw new IllegalArgumentException("${fieldName} contains duplicate barcodes")
        }

        normalized.each { barcode ->
            if( barcode.size() != expectedLength ) {
                throw new IllegalArgumentException(
                    "${fieldName} contains barcode '${barcode}' with length ${barcode.size()}; expected ${expectedLength} nt"
                )
            }
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

    private static File writeSbGroupMap(final File derivedDir, final List<Map> samples, final String modality) {
        final File file = new File(derivedDir, "${modality}_sb_group_map.tsv")
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

    private static String parseDnaTagmentation(final Map dnaConfig, final String sampleId) {
        final String mode = (dnaConfig.tagmentation ?: TAGMENTATION_SINGLE).toString().trim().toLowerCase()
        if( !(mode in DNA_TAGMENTATION_MODES) ) {
            throw new IllegalArgumentException(
                "samples.${sampleId}.dna.tagmentation must be one of: single, dual"
            )
        }
        return mode
    }

    private static int dnaSbBarcodeLength(final String tagmentation) {
        return tagmentation == TAGMENTATION_DUAL ? DNA_DUAL_SB_BARCODE_LENGTH : DNA_SINGLE_SB_BARCODE_LENGTH
    }

    private static Map dnaTagDefaultsForTagmentation(final Map defaults, final String tagmentation) {
        if( tagmentation == TAGMENTATION_SINGLE ) {
            return [
                sample  : [
                    bc_len            : defaults.dna.sample.bc_len as int,
                    bc_start          : defaults.dna.sample.bc_start as int,
                    hd                : defaults.dna.sample.hd as int,
                    reverse_complement: toDirection(defaults.dna.sample.reverse_complement, 'barcode_defaults.dna.sample.reverse_complement'),
                    index_read        : 'i2',
                ],
                modality: [
                    bc_len            : defaults.dna.modality.bc_len as int,
                    bc_start          : defaults.dna.modality.bc_start as int,
                    hd                : defaults.dna.modality.hd as int,
                    reverse_complement: toDirection(defaults.dna.modality.reverse_complement, 'barcode_defaults.dna.modality.reverse_complement'),
                    index_read        : 'i2',
                ],
            ]
        }

        return [
            sample  : [
                bc_len            : 3,
                bc_start          : 0,
                hd                : 0,
                reverse_complement: 'fw',
                index_read        : 'i1',
            ],
            modality: [
                bc_len            : 8,
                bc_start          : 3,
                hd                : 1,
                reverse_complement: 'fw',
                index_read        : 'i1',
            ],
        ]
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

    private static Map resolveReads(
        final File baseDir,
        final Map modalityConfig,
        final String sampleId,
        final String modality,
        final List<String> readNames
    ) {
        final Map reads = asMap(modalityConfig.reads, "samples.${sampleId}.${modality}.reads")
        return readNames.collectEntries { readName ->
            final String fieldName = "samples.${sampleId}.${modality}.reads.${readName}"
            [(readName): resolveExistingPath(baseDir, requireString(reads[readName], fieldName))]
        }
    }

    private static Map resolveRuntime(final Map parsed, final File baseDir) {
        final Map runtime = asMap(parsed.runtime, 'runtime')
        return [
            env_prefix: resolvePath(
                baseDir,
                requireString(runtime.env_prefix, 'runtime.env_prefix')
            ),
            tmpdir: resolvePath(
                baseDir,
                requireString(runtime.tmpdir, 'runtime.tmpdir')
            ),
        ]
    }

    private static Map resolveReferences(final Map parsed, final File baseDir) {
        final Map references = asMap(parsed.references, 'references')
        final String root = resolvePath(
            baseDir,
            requireString(references.root, 'references.root')
        )
        final File rootDir = new File(root)
        final String rnaRefDir = resolveOptionalPath(baseDir, references.rna_ref_dir)

        return [
            species                    : requireString(references.species, 'references.species').toLowerCase(),
            root                       : rootDir.canonicalPath,
            ligation_barcode_whitelist : resolvePath(
                baseDir,
                requireString(references.ligation_barcode_whitelist, 'references.ligation_barcode_whitelist')
            ),
            rna_ref_dir                : rnaRefDir,
            rna_chrom_sizes            : rnaRefDir ? new File(rnaRefDir, 'chrNameLength.txt').canonicalPath : null,
            dna_ref_dir                : resolveOptionalPath(baseDir, references.dna_ref_dir),
            dna_bwa_reference          : null,
            dna_blacklist_bed          : resolveOptionalPath(baseDir, references.dna_blacklist_bed),
            dna_chrom_sizes            : resolveOptionalPath(baseDir, references.dna_chrom_sizes),
            dna_effective_genome_size  : references.dna_effective_genome_size?.toString()?.trim(),
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

    private static String resolvePath(final File baseDir, final String rawPath) {
        final File resolved = new File(rawPath).isAbsolute() ? new File(rawPath) : new File(baseDir, rawPath)
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
