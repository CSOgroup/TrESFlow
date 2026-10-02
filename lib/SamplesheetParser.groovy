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
    private static final List<String> FASTQ_SUFFIXES = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
    private static final java.util.regex.Pattern CONTROL_CHARACTER_PATTERN =
        java.util.regex.Pattern.compile('[\\x00-\\x1F\\x7F]')

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
        final Map runtime = resolveRuntime(parsed, baseDir, options)
        final Map references = resolveReferences(parsed, baseDir)
        final Map resolved = parseUnified(parsed, baseDir, libraryName, options, references)
        final List<Map> samples = resolved.samples
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
            sb_identity_warnings: resolved.warnings,
        ]
    }

    static String normalizeOligoIndex(final Object value, final String fieldName) {
        if( value instanceof Boolean || !(value instanceof String || value instanceof Byte || value instanceof Short || value instanceof Integer || value instanceof Long || value instanceof BigInteger) ) {
            throw new IllegalArgumentException("${fieldName} requires positive integers or decimal strings; received '${value}'")
        }
        final String text = value.toString()
        if( !(text ==~ /[0-9]+/) || new BigInteger(text) <= 0 ) {
            throw new IllegalArgumentException("${fieldName} requires positive decimal indices; received '${value}'")
        }
        return new BigInteger(text).toString().padLeft(2, '0')
    }

    static Map loadOligoLookup(final File file) {
        if( !file.isFile() ) throw new IllegalArgumentException("SB oligo lookup not found: ${file}")
        final List<String> lines = file.readLines()
        if( !lines || lines[0] != 'oligo_index\trna_and_single_dna_sb\tdual_dna_sb' ) {
            throw new IllegalArgumentException("Invalid SB oligo lookup header: ${file}")
        }
        final Map rows = [:]
        final Map four = [:]
        final Map dual = [:]
        lines.drop(1).eachWithIndex { line, n ->
            final String field = "SB oligo lookup ${file}:${n + 2}"
            final List parts = line.split('\t', -1).toList()
            if( parts.size() != 3 ) throw new IllegalArgumentException("${field}: expected three nonempty columns")
            final String index = normalizeOligoIndex(parts[0], field)
            if( rows.containsKey(index) ) throw new IllegalArgumentException("${field}: duplicate normalized index '${index}'")
            if( !(parts[1] ==~ /[ACGT]{4}/) || !(parts[2] == '-' || parts[2] ==~ /[ACGT]{3}/) ) {
                throw new IllegalArgumentException("${field}: requires uppercase A/C/G/T four-base and three-base sequences (or '-' for unavailable dual DNA)")
            }
            if( four.containsKey(parts[1]) || (parts[2] != '-' && dual.containsKey(parts[2])) ) {
                throw new IllegalArgumentException("${field}: ambiguous duplicate chemistry sequence")
            }
            rows[index] = [parts[1], parts[2]]
            four[parts[1]] = index
            if( parts[2] != '-' ) dual[parts[2]] = index
        }
        if( rows.isEmpty() ) throw new IllegalArgumentException("SB oligo lookup has no rows: ${file}")
        return [rows: rows, four: four, dual: dual]
    }

    private static List<String> requireOligoIndices(final Object value, final String field, final Map lookup) {
        if( !(value instanceof List) || value.isEmpty() ) throw new IllegalArgumentException("${field} must be a non-empty list")
        final List<String> indices = value.collect { normalizeOligoIndex(it, field) }
        if( indices.toSet().size() != indices.size() ) throw new IllegalArgumentException("${field} contains duplicate normalized indices")
        indices.each { if( !lookup.rows.containsKey(it) ) throw new IllegalArgumentException("${field}: unknown oligo index '${it}'") }
        return indices
    }

    private static String requireIdentityName(final Object value, final String field) {
        final String name = requireString(value, field)
        if( name != value.toString() || !(name ==~ /[A-Za-z0-9_.-]+/) || name in ['.', '..', 'Unknown'] ) {
            throw new IllegalArgumentException("${field}: unsupported identity/output name '${name}' (use letters, digits, underscores, dots or hyphens)")
        }
        return name
    }

    private static void registerSerializedName(final Map seen, final String name, final List owner, final String kind) {
        if( seen.containsKey(name) && seen[name] != owner ) {
            throw new IllegalArgumentException("Serialized ${kind} collision '${name}': ${seen[name]} vs ${owner}")
        }
        seen[name] = owner
    }

    private static String injectedSbBase(final String firstPass) {
        if( firstPass == 'first_pass' ) return ''
        if( firstPass ==~ /first_pass_withBC_[ACGT]/ ) return firstPass[-1]
        throw new IllegalArgumentException("Unsupported SB upstream tag format '${firstPass}'")
    }

    private static Map parseUnified(
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

        final File lookupFile = new File((options.sb_oligo_lookup ?: 'assets/sb_oligo_lookup.v1.tsv').toString())
        final Map lookup = loadOligoLookup(lookupFile)
        final String lookupDigest = java.security.MessageDigest.getInstance('SHA-256')
            .digest(lookupFile.bytes).encodeHex().toString()
        final Map namespaces = [:]
        final Map outputs = [:]
        final List<Map> samples = []
        final List<String> warnings = []

        ((Map) parsed.samples).each { rawSampleId, rawSampleConfig ->
            final String sampleId = requireIdentityName(rawSampleId, 'samples.<sample_id>')
            final Map sampleConfig = asMap(rawSampleConfig, "samples.${sampleId}")
            final Map groupsConfig = asMap(sampleConfig.groups, "samples.${sampleId}.groups")
            if( groupsConfig.isEmpty() ) {
                throw new IllegalArgumentException("samples.${sampleId}.groups must not be empty")
            }

            final Map rnaConfig = sampleConfig.rna ? asMap(sampleConfig.rna, "samples.${sampleId}.rna") : null
            final Map dnaConfig = sampleConfig.dna ? asMap(sampleConfig.dna, "samples.${sampleId}.dna") : null
            final boolean hasRna = rnaConfig != null
            final boolean hasDna = dnaConfig != null

            if( hasDna && dnaConfig.containsKey('mark_barcodes') ) {
                throw new IllegalArgumentException(
                    "samples.${sampleId}.dna.mark_barcodes is no longer supported; " +
                    "define mark_barcodes under each DNA group at " +
                    "samples.${sampleId}.groups.<group>.mark_barcodes"
                )
            }

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
                dnaTagmentation,
                lookup
            )

            groupsConfig.keySet().each { group ->
                registerSerializedName(namespaces, "${sampleId}_${group}", [sampleId, group.toString()], 'cell namespace')
            }
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
                    normalizedGroups.dna_mark_barcodes as LinkedHashMap<String, Map<String, String>>,
                    normalizedGroups.dna_source_summary as String,
                    dnaTagmentation,
                    dnaConfig,
                    defaults,
                    ligationWhitelist,
                    references
                )
            }
            warnings.addAll(normalizedGroups.warnings)
            samples.findAll { it.id == sampleId }.each { row ->
                row.group_sources = normalizedGroups[row.modality + '_sources']
                row.group_identities = normalizedGroups[row.modality + '_identities']
            }
        }

        samples.each { row ->
            row.cell_id_version = 'full-cell-v2'
            row.sb_lookup_sha256 = lookupDigest
            row.sb_injected_base = injectedSbBase(row.sample_first_pass as String)
            row.cell_identity_records = row.group_definitions.collectMany { group, sequences ->
                sequences.collect { sequence ->
                    [sample: row.id, group: group, modality: row.modality,
                     chemistry: row.modality == MODALITY_RNA ? 'rna' : "dna_${row.dna_tagmentation}",
                     oligo_index: lookup[row.modality == MODALITY_DNA && row.dna_tagmentation == TAGMENTATION_DUAL ? 'dual' : 'four'][sequence],
                     sb_bc: sequence, input_source: row.group_identities[group][sequence].input_source,
                     sb_index: row.group_identities[group][sequence].sb_index]
                }
            }
            row.sb_mapping_sha256 = java.security.MessageDigest.getInstance('SHA-256')
                .digest(groovy.json.JsonOutput.toJson(row.cell_identity_records).getBytes('UTF-8')).encodeHex().toString()
            row.split_targets = [:]
            row.group_definitions.keySet().each { group ->
                final List marks = row.modality == MODALITY_DNA ? row.mark_barcodes[group].keySet().toList() : ['RNA']
                marks.each { mark ->
                    final String stem = row.modality == MODALITY_DNA ? "${row.id}_${group}_${mark}" : "${row.id}_${group}"
                    registerSerializedName(outputs, "${row.modality}:${stem}", [row.id, group, mark], 'output name')
                    row.split_targets[stem] = [sample: row.id, group: group, mark: mark, sample_group: "${row.id}_${group}"]
                }
            }
            row.remove('group_sources')
            row.remove('group_identities')
        }
        writeDerivedText(new File(derivedDir, 'sb_oligo_lookup.v1.tsv'), lookupFile.getText('UTF-8'))
        writeDerivedText(new File(derivedDir, 'cell_identity_version.txt'), "full-cell-v2\nCB=XI=<sample>_<group>_<sb_index>_<L1L2L3>\nsb_index=physical oligo_index unless explicitly paired; SB=actual corrected chemistry sequence\nlookup_sha256=${lookupDigest}\n")
        final List<Map> attached = attachDerivedArtifacts(derivedDir, samples)
        final File infoDir = options.outdir?.toString()?.trim() ? derivedDir.parentFile : derivedDir
        writeDerivedText(new File(infoDir, 'sb_identity_warnings.txt'), warnings
            ? warnings.collect { delimitIdentityWarning(it) }.join('\n\n') + '\n'
            : 'No SB identity warnings.\n')
        final List<String> audit = ['sample\tsb_group\tsb_bc\toligo_index\tmodality\tchemistry\tinput_source\tsb_injected_base\tsb_index']
        samples.each { row ->
            row.cell_identity_records.each { r ->
                audit << "${r.sample}\t${r.group}\t${r.sb_bc}\t${r.oligo_index}\t${r.modality}\t${r.chemistry}\t${r.input_source}\t${row.sb_injected_base ?: '-'}\t${r.sb_index}"
            }
        }
        writeDerivedText(new File(infoDir, 'sb_physical_to_logical.tsv'), audit.join('\n') + '\n')
        return [samples: attached, warnings: warnings]
    }

    static String delimitIdentityWarning(final String warning) {
        return "================ SB IDENTITY WARNING ================\n${warning}\n====================================================="
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
        final LinkedHashMap<String, Map<String, String>> markBarcodesByGroup,
        final String sbSourceSummary,
        final String tagmentation,
        final Map dnaConfig,
        final Map defaults,
        final String ligationWhitelist,
        final Map references
    ) {
        final Map reads = resolveDnaReads(baseDir, dnaConfig, sampleId, tagmentation)
        final Map dnaTagDefaults = dnaTagDefaultsForTagmentation(defaults, tagmentation)

        return [
            id                          : sampleId,
            modality                    : MODALITY_DNA,
            dna_tagmentation            : tagmentation,
            i1                          : reads.i1,
            i2                          : reads.i2,
            i2_implicit                 : reads.i2_implicit,
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
            mark_barcodes               : markBarcodesByGroup,
        ]
    }

    private static List<Map> attachDerivedArtifacts(final File derivedDir, final List<Map> samples) {
        final List<Map> rnaRows = samples.findAll { it.modality == MODALITY_RNA }
        final List<Map> dnaRows = samples.findAll { it.modality == MODALITY_DNA }
        final File rnaSbGroupMapFile = rnaRows ? writeSbGroupMap(derivedDir, rnaRows, MODALITY_RNA) : null
        final File dnaSbGroupMapFile = dnaRows ? writeSbGroupMap(derivedDir, dnaRows, MODALITY_DNA) : null
        final File dnaMoMapFile = dnaRows ? writeDnaMoMap(derivedDir, dnaRows) : null
        final Map<String, String> dnaWhitelistPaths = dnaRows ? writeDnaModalityWhitelists(derivedDir, dnaRows) : [:]
        writeInputFastqProvenance(derivedDir, samples)

        samples.each { row ->
            row.sb_group_map = row.modality == MODALITY_RNA
                ? rnaSbGroupMapFile.canonicalPath
                : dnaSbGroupMapFile.canonicalPath
            // Preserve the samplesheet's explicit group identities as report
            // metadata before removing the internal barcode-definition map.
            row.samplesheet_groups = new ArrayList<String>(
                (row.group_definitions as Map).keySet().collect { it.toString() }
            )
            row.remove('group_definitions')
            if( row.modality == MODALITY_DNA ) {
                row.mo_map = dnaMoMapFile.canonicalPath
                row.modality_whitelist = dnaWhitelistPaths[row.id]
                row.remove('mark_barcodes')
            }
        }

        return samples
    }

    private static File writeInputFastqProvenance(final File derivedDir, final List<Map> samples) {
        final File file = new File(derivedDir, 'input_fastq_provenance.tsv')
        final List<String> lines = [
            'sample\tmodality\tread_set_index\tread_role\tcanonical_input_path\tinput_origin'
        ]

        samples.each { row ->
            final List<String> roles = row.modality == MODALITY_DNA
                ? ['i1', 'i2', 'r1', 'r2']
                : ['i1', 'r1', 'r2']
            roles.each { role ->
                (row[role] as List<String>).eachWithIndex { path, index ->
                    final String origin = row.modality == MODALITY_DNA && role == 'i2' && row.i2_implicit
                        ? 'implicit_i1_fallback'
                        : 'explicit'
                    lines << "${row.id}\t${row.modality}\t${index + 1}\t${role}\t${path}\t${origin}"
                }
            }
        }

        writeDerivedText(file, lines.join('\n') + '\n')
        return file
    }

    private static Map parseGroups(
        final Map groupsConfig,
        final String sampleId,
        final boolean hasRna,
        final boolean hasDna,
        final String dnaTagmentation,
        final Map lookup
    ) {
        final LinkedHashMap<String, List<String>> rnaGroups = new LinkedHashMap<>()
        final LinkedHashMap<String, List<String>> dnaGroups = new LinkedHashMap<>()
        final LinkedHashMap<String, Map<String, String>> dnaMarkBarcodes = new LinkedHashMap<>()
        final LinkedHashMap<String, String> rnaSources = new LinkedHashMap<>()
        final LinkedHashMap<String, String> dnaSources = new LinkedHashMap<>()
        final Map identities = [rna: [:], dna: [:]]
        final List<String> warnings = []
        final List<String> selectors = ['sb_oligo_indices', 'rna_sb_oligo_indices', 'dna_sb_oligo_indices',
                                        'sb_barcodes', 'rna_sb_barcodes', 'dna_sb_barcodes']

        groupsConfig.each { rawGroupName, rawGroupConfig ->
            final String groupName = requireIdentityName(rawGroupName, "samples.${sampleId}.groups.<group>")
            final String field = "samples.${sampleId}.groups.${groupName}"
            final Map groupConfig = asMap(rawGroupConfig, field)
            final boolean paired = groupConfig.containsKey('sb_oligo_pairings')
            final boolean shared = groupConfig.containsKey('sb_oligo_indices')
            if( paired && selectors.any { groupConfig.containsKey(it) } ) {
                throw new IllegalArgumentException("${field}: sb_oligo_pairings cannot be combined with other barcode-selection fields")
            }
            if( shared && selectors.findAll { it != 'sb_oligo_indices' }.any { groupConfig.containsKey(it) } ) {
                throw new IllegalArgumentException("${field}: sb_oligo_indices cannot be combined with other barcode-selection fields")
            }
            final List<String> sharedIndices = shared ? requireOligoIndices(groupConfig.sb_oligo_indices, "${field}.sb_oligo_indices", lookup) : []
            final Map selected = paired
                ? resolvePairings(groupConfig.sb_oligo_pairings, field, hasRna, hasDna, dnaTagmentation, lookup)
                : [rna: [:], dna: [:]]
            if( !paired ) {
                [rna: hasRna, dna: hasDna].each { modality, present ->
                    final String indexField = "${modality}_sb_oligo_indices"
                    final String sequenceField = "${modality}_sb_barcodes"
                    final boolean specificIndices = groupConfig.containsKey(indexField)
                    if( specificIndices && (groupConfig.containsKey(sequenceField) ||
                        (groupConfig.containsKey('sb_barcodes') && (modality == MODALITY_RNA || dnaTagmentation == TAGMENTATION_SINGLE))) ) {
                        throw new IllegalArgumentException("${field}: conflicting selectors ${indexField} and sequence-input fields")
                    }
                    // Validate explicit index lists even when a modality is absent.
                    final List<String> indices = specificIndices ? requireOligoIndices(groupConfig[indexField], "${field}.${indexField}", lookup) : sharedIndices
                    final boolean hasSequences = groupConfig.containsKey(sequenceField) ||
                        (groupConfig.containsKey('sb_barcodes') && (modality == MODALITY_RNA || dnaTagmentation == TAGMENTATION_SINGLE))
                    if( present && (shared || specificIndices || hasSequences) ) {
                        final boolean dual = modality == MODALITY_DNA && dnaTagmentation == TAGMENTATION_DUAL
                        final String source = "${field}.${shared ? 'sb_oligo_indices' : indexField}"
                        final Map selection = shared || specificIndices
                            ? [value: indices.collect { physicalSequence(it, dual, source, lookup) }, fieldName: source]
                            : (modality == MODALITY_RNA ? selectRnaSbBarcodes(groupConfig, sampleId, groupName)
                                : selectDnaSbBarcodes(groupConfig, sampleId, groupName, dnaTagmentation))
                        requireBarcodeList(selection.value, selection.fieldName, dual ? 3 : 4).each { sequence ->
                            final String index = physicalIndex(sequence, dual, selection.fieldName, lookup)
                            selected[modality][sequence] = [sb_index: index, input_source: selection.fieldName]
                        }
                    }
                }
            }

            if( !selected.rna.isEmpty() ) {
                rnaGroups[groupName] = selected.rna.keySet().toList()
                rnaSources[groupName] = selected.rna.values().first().input_source
                identities.rna[groupName] = selected.rna
            }
            final boolean hasDnaBarcodes = !selected.dna.isEmpty()
            final boolean hasMarks = groupConfig.containsKey('mark_barcodes')
            if( hasDna && hasMarks && !hasDnaBarcodes ) {
                if( !paired && dnaTagmentation == TAGMENTATION_DUAL ) {
                    selectDnaSbBarcodes(groupConfig, sampleId, groupName, dnaTagmentation)
                }
                throw new IllegalArgumentException("${field}.mark_barcodes requires a DNA sample-barcode field (dna_sb_barcodes, dna_sb_oligo_indices or a DNA pairing)")
            }
            if( hasDnaBarcodes && !hasMarks ) {
                throw new IllegalArgumentException("Missing required field: ${field}.mark_barcodes for DNA group")
            }
            if( hasDnaBarcodes ) {
                dnaGroups[groupName] = selected.dna.keySet().toList()
                dnaSources[groupName] = selected.dna.values().first().input_source
                identities.dna[groupName] = selected.dna
                dnaMarkBarcodes[groupName] = parseMarkBarcodes(groupConfig.mark_barcodes, sampleId, groupName)
            }

            if( paired ) {
                final Set labels = (selected.rna.values() + selected.dna.values()).collect { it.sb_index }.toSet()
                labels.each { label ->
                    final Map physical = [:]
                    ['rna', 'dna'].each { modality ->
                        final def entry = selected[modality].find { sequence, record -> record.sb_index == label }
                        if( entry ) physical[modality] = [index: physicalIndex(entry.key, modality == 'dna' && dnaTagmentation == TAGMENTATION_DUAL, field, lookup), sequence: entry.key]
                    }
                    if( physical.values().any { it.index != label } ) {
                        final String rna = physical.rna ? "${physical.rna.index} (${physical.rna.sequence})" : 'omitted'
                        final String dna = physical.dna ? "${physical.dna.index} (${physical.dna.sequence})" : 'omitted'
                        warnings << "Sample '${sampleId}', group '${groupName}': actual RNA oligo ${rna}; actual DNA oligo ${dna}; DNA chemistry ${dnaTagmentation ?: 'not selected'}; shared sb_index ${label}. " +
                            (physical.rna && physical.dna
                                ? (physical.rna.index != physical.dna.index
                                    ? 'The samplesheet explicitly declares these different physical oligos to represent one biological partition.'
                                    : 'The samplesheet explicitly declares these physical oligos to represent one biological partition under a remapped identity label.')
                                : 'The samplesheet explicitly assigns this physical oligo to the biological partition identified by sb_index; the other modality is omitted.') +
                            ' Pairing changes CB/XI identity labels and retains the actual chemistry-specific SB sequences.'
                    }
                }
            } else if( !selected.rna.isEmpty() && !selected.dna.isEmpty() ) {
                final Set rnaIndices = selected.rna.values().collect { it.sb_index }.toSet()
                final Set dnaIndices = selected.dna.values().collect { it.sb_index }.toSet()
                if( rnaIndices != dnaIndices ) {
                    warnings << "Sample '${sampleId}', group '${groupName}': independent RNA oligos ${describeSelection(selected.rna, false, lookup)} and DNA oligos ${describeSelection(selected.dna, dnaTagmentation == TAGMENTATION_DUAL, lookup)} differ; DNA chemistry ${dnaTagmentation}. " +
                        'No explicit pairing was declared. Each modality retains its actual oligo index as sb_index, so unmatched full identifiers will differ. Correspondence is never inferred from list order.'
                }
            }
        }

        [rna: rnaGroups, dna: dnaGroups].each { modality, groups ->
            if( (modality == MODALITY_RNA ? hasRna : hasDna) && groups.isEmpty() ) {
                throw new IllegalArgumentException("samples.${sampleId} has ${modality.toUpperCase()} reads but no group with ${modality}_sb_barcodes${modality == MODALITY_DNA ? ' and mark_barcodes' : ''}")
            }
            validateNoGroupBarcodeCollisions(sampleId, modality.toUpperCase(), groups)
        }
        // Same physical oligo in different modalities has unambiguous routing;
        // group namespaces remain explicit and distinct.
        rnaGroups.each { rnaGroup, rnaSequences ->
            dnaGroups.each { dnaGroup, dnaSequences ->
                if( rnaGroup != dnaGroup ) {
                    final Set overlap = rnaSequences.collect { lookup.four[it] }.toSet().intersect(
                        dnaSequences.collect { (dnaTagmentation == TAGMENTATION_DUAL ? lookup.dual : lookup.four)[it] }.toSet())
                    if( overlap ) warnings << "Sample '${sampleId}': physical oligo indices ${overlap} occur in RNA group '${rnaGroup}' and DNA group '${dnaGroup}' (DNA chemistry ${dnaTagmentation}). Routing is unambiguous within each modality; their full identifiers differ because their group namespaces differ."
                }
            }
        }
        return [rna_sources: rnaSources, dna_sources: dnaSources, rna: rnaGroups, dna: dnaGroups,
                rna_identities: identities.rna, dna_identities: identities.dna, warnings: warnings,
                dna_mark_barcodes: dnaMarkBarcodes, rna_source_summary: summarizeSources(rnaSources.values()),
                dna_source_summary: summarizeSources(dnaSources.values())]
    }

    private static String physicalSequence(final String index, final boolean dual, final String field, final Map lookup) {
        if( !lookup.rows.containsKey(index) ) throw new IllegalArgumentException("${field}: unknown oligo index '${index}'")
        final String sequence = lookup.rows[index][dual ? 1 : 0]
        if( sequence == '-' ) throw new IllegalArgumentException("${field}: physical oligo '${index}' has no available dual-DNA sequence. Select an available physical dna_oligo_index/dna_sb_oligo_indices or use RNA/single DNA; sb_index is only an identity label. No dual sequence will be inferred or substituted.")
        return sequence
    }

    private static String physicalIndex(final String sequence, final boolean dual, final String field, final Map lookup) {
        final String index = (dual ? lookup.dual : lookup.four)[sequence]
        if( !index ) throw new IllegalArgumentException("${field}: Unknown SB sequence '${sequence}' for ${dual ? 'dual DNA' : 'RNA/single DNA'} in oligo lookup")
        return index
    }

    private static String describeSelection(final Map selected, final boolean dual, final Map lookup) {
        return selected.keySet().collect { "${physicalIndex(it, dual, 'SB selection', lookup)} (${it})" }.join(', ')
    }

    private static Map resolvePairings(final Object value, final String field, final boolean hasRna,
                                      final boolean hasDna, final String dnaTagmentation, final Map lookup) {
        final String prefix = "${field}.sb_oligo_pairings"
        if( !(value instanceof List) || value.isEmpty() ) throw new IllegalArgumentException("${prefix} must be a non-empty list")
        final Map selected = [rna: [:], dna: [:]]
        final Map logical = [rna: [:], dna: [:]]
        value.eachWithIndex { raw, n ->
            final String entryField = "${prefix}[${n + 1}]"
            final Map entry = asMap(raw, entryField)
            final List allowed = ['sb_index', 'rna_oligo_index', 'dna_oligo_index', 'rna_sb_barcode', 'dna_sb_barcode']
            if( entry.keySet().any { !(it in allowed) } ) throw new IllegalArgumentException("${entryField}: unknown pairing fields ${entry.keySet().findAll { !(it in allowed) }}; use sb_index and physical modality selectors")
            final String label = normalizeOligoIndex(entry.sb_index, "${entryField}.sb_index")
            boolean anyModality = false
            [rna: hasRna, dna: hasDna].each { modality, present ->
                final String indexField = "${modality}_oligo_index"
                final String sequenceField = "${modality}_sb_barcode"
                final List fields = [indexField, sequenceField].findAll { entry.containsKey(it) }
                if( fields.size() > 1 ) throw new IllegalArgumentException("${entryField}: exactly one ${modality} selector is required; ${indexField} conflicts with ${sequenceField}")
                if( fields ) {
                    anyModality = true
                    if( !present ) throw new IllegalArgumentException("${entryField}: ${modality} selector requires a sample ${modality} reads block")
                    final boolean dual = modality == MODALITY_DNA && dnaTagmentation == TAGMENTATION_DUAL
                    final String source = "${entryField}.${fields[0]}"
                    final String sequence = fields[0] == indexField
                        ? physicalSequence(normalizeOligoIndex(entry[indexField], source), dual, source, lookup)
                        : requireBarcodeList([entry[sequenceField]], source, dual ? 3 : 4).first()
                    physicalIndex(sequence, dual, source, lookup)
                    if( selected[modality].containsKey(sequence) ) throw new IllegalArgumentException("${entryField}: physical ${modality} barcode '${sequence}' is reused for multiple partitions")
                    if( logical[modality].containsKey(label) ) throw new IllegalArgumentException("${entryField}: sb_index '${label}' allows at most one physical oligo per modality (${modality})")
                    selected[modality][sequence] = [sb_index: label, input_source: source]
                    logical[modality][label] = sequence
                }
            }
            if( !anyModality ) throw new IllegalArgumentException("${entryField}: at least one modality selector is required")
        }
        return selected
    }

    private static LinkedHashMap<String, String> parseMarkBarcodes(
        final Object value,
        final String sampleId,
        final String groupName
    ) {
        final String fieldPrefix = "samples.${sampleId}.groups.${groupName}.mark_barcodes"
        final Map marksConfig = asMap(value, fieldPrefix)
        if( marksConfig.isEmpty() ) {
            throw new IllegalArgumentException("${fieldPrefix} must not be empty")
        }

        final LinkedHashMap<String, String> normalized = new LinkedHashMap<>()
        final Map<String, String> barcodeToMark = [:]

        marksConfig.each { rawMarkName, rawBarcode ->
            final String markName = requireIdentityName(rawMarkName, "${fieldPrefix}.<mark>")
            final String barcode = requireString(rawBarcode, "${fieldPrefix}.${markName}")
            if( barcodeToMark.containsKey(barcode) && barcodeToMark[barcode] != markName ) {
                throw new IllegalArgumentException(
                    "Duplicate DNA modality barcode '${barcode}' for sample '${sampleId}' group '${groupName}': " +
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

        root.mkdirs()
        return root
    }

    private static void writeDerivedText(final File file, final String text) {
        // Keep unchanged staged contracts stable for Nextflow's default cache.
        if( !file.isFile() || file.getText('UTF-8') != text ) file.setText(text, 'UTF-8')
    }

    private static File writeSbGroupMap(final File derivedDir, final List<Map> samples, final String modality) {
        final File file = new File(derivedDir, "${modality}_sb_group_map.tsv")
        final List<String> lines = ['sample\tsb_group\tsb_bc\toligo_index\tmodality\tchemistry\tinput_source\tsb_injected_base\tsb_index']

        samples.collect { it.id }.unique().each { sampleId ->
            final Map row = samples.find { it.id == sampleId }
            row.group_definitions.each { groupName, barcodes ->
                barcodes.each { barcode ->
                    final Map record = row.cell_identity_records.find { it.group == groupName && it.sb_bc == barcode }
                    lines << "${sampleId}\t${groupName}\t${barcode}\t${record.oligo_index}\t${modality}\t${record.chemistry}\t${record.input_source}\t${row.sb_injected_base ?: '-'}\t${record.sb_index}"
                }
            }
        }

        writeDerivedText(file, lines.join('\n') + '\n')
        return file
    }

    private static String parseDnaTagmentation(final Map dnaConfig, final String sampleId) {
        final String mode = requireString(
            dnaConfig.tagmentation,
            "samples.${sampleId}.dna.tagmentation"
        ).toLowerCase()
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
                final Map groupMarkBarcodes = (row.mark_barcodes as Map)[groupName] as Map
                groupMarkBarcodes.each { markName, barcode ->
                    lines << "${row.id}\t${groupName}\t${markName}\t${barcode}"
                }
            }
        }

        writeDerivedText(file, lines.join('\n') + '\n')
        return file
    }

    private static Map<String, String> writeDnaModalityWhitelists(final File derivedDir, final List<Map> dnaRows) {
        final File whitelistDir = new File(derivedDir, 'dna_modality_whitelists')
        whitelistDir.mkdirs()

        final Map<String, String> out = [:]
        dnaRows.each { row ->
            final File file = new File(whitelistDir, "${row.id}.txt")
            final LinkedHashSet<String> modalityBarcodes = new LinkedHashSet<>()
            (row.mark_barcodes as Map).values().each { groupMarks ->
                (groupMarks as Map).values().each { barcode -> modalityBarcodes.add(barcode.toString()) }
            }
            writeDerivedText(file, modalityBarcodes.join('\n') + '\n')
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
        final Map resolved = readNames.collectEntries { readName ->
            final String fieldName = "samples.${sampleId}.${modality}.reads.${readName}"
            [(readName): resolveReadPaths(baseDir, reads[readName], fieldName)]
        }
        validateReadSetContract(sampleId, modality, resolved, false)
        return resolved
    }

    private static Map resolveDnaReads(
        final File baseDir,
        final Map dnaConfig,
        final String sampleId,
        final String tagmentation
    ) {
        final List<String> requiredReadNames = tagmentation == TAGMENTATION_DUAL
            ? ['i1', 'r1', 'r2']
            : ['i1', 'i2', 'r1', 'r2']
        final Map reads = asMap(dnaConfig.reads, "samples.${sampleId}.dna.reads")
        final Map resolved = requiredReadNames.collectEntries { readName ->
            final String fieldName = "samples.${sampleId}.dna.reads.${readName}"
            [(readName): resolveReadPaths(baseDir, reads[readName], fieldName)]
        }

        if( tagmentation == TAGMENTATION_DUAL ) {
            final boolean i2Explicit = reads.containsKey('i2')
            resolved.i2 = i2Explicit
                ? resolveReadPaths(baseDir, reads.i2, "samples.${sampleId}.dna.reads.i2")
                : new ArrayList<String>(resolved.i1 as List<String>)
            resolved.i2_implicit = !i2Explicit
        }
        else {
            resolved.i2_implicit = false
        }

        validateReadSetContract(sampleId, MODALITY_DNA, resolved, resolved.i2_implicit as boolean)

        return resolved
    }

    private static List<String> resolveReadPaths(
        final File baseDir,
        final Object rawValue,
        final String fieldName
    ) {
        final List rawEntries
        if( rawValue instanceof List ) {
            if( rawValue.isEmpty() ) {
                throw new IllegalArgumentException("${fieldName} must contain at least one FASTQ path")
            }
            rawEntries = rawValue as List
        }
        else {
            if( rawValue == null ) {
                throw new IllegalArgumentException("Missing required field: ${fieldName}")
            }
            if( !(rawValue instanceof CharSequence) ) {
                throw new IllegalArgumentException(
                    "${fieldName} must be a FASTQ path string or a list of FASTQ path strings"
                )
            }
            final String scalar = rawValue.toString()
            if( !scalar.trim() ) {
                throw new IllegalArgumentException("${fieldName} must contain at least one FASTQ path")
            }
            rawEntries = Arrays.asList(scalar.split(',', -1))
        }

        final List<String> resolved = []
        rawEntries.eachWithIndex { rawEntry, index ->
            if( rawEntry == null || rawEntry instanceof Map || rawEntry instanceof List ) {
                throw new IllegalArgumentException(
                    "${fieldName} entry ${index + 1} must be a FASTQ path string"
                )
            }
            final String rawPath = rawEntry.toString()
            if( CONTROL_CHARACTER_PATTERN.matcher(rawPath).find() ) {
                throw new IllegalArgumentException(
                    "${fieldName} entry ${index + 1} contains a control character"
                )
            }
            final String path = rawPath.trim()
            if( !path ) {
                throw new IllegalArgumentException(
                    "${fieldName} contains an empty FASTQ path at entry ${index + 1}"
                )
            }
            final File candidate = new File(path).isAbsolute() ? new File(path) : new File(baseDir, path)
            if( !candidate.exists() ) {
                throw new IllegalArgumentException(
                    "${fieldName} entry ${index + 1} FASTQ not found: ${candidate}"
                )
            }
            if( !candidate.isFile() ) {
                throw new IllegalArgumentException(
                    "${fieldName} entry ${index + 1} is not a regular file: ${candidate}"
                )
            }
            final String canonicalPath = candidate.canonicalPath
            final String lowercaseName = candidate.name.toLowerCase()
            if( !FASTQ_SUFFIXES.any { suffix -> lowercaseName.endsWith(suffix) } ) {
                throw new IllegalArgumentException(
                    "${fieldName} entry ${index + 1} has an invalid FASTQ suffix: ${candidate}; " +
                    "expected .fastq, .fq, .fastq.gz, or .fq.gz"
                )
            }
            resolved << canonicalPath
        }

        if( resolved.toSet().size() != resolved.size() ) {
            throw new IllegalArgumentException(
                "${fieldName} contains duplicate canonical FASTQ paths"
            )
        }
        return resolved
    }

    private static void validateReadSetContract(
        final String sampleId,
        final String modality,
        final Map reads,
        final boolean implicitDualI2
    ) {
        final List<String> roles = modality == MODALITY_DNA
            ? ['i1', 'i2', 'r1', 'r2']
            : ['i1', 'r1', 'r2']
        final Map<String, Integer> counts = roles.collectEntries { role ->
            [(role): (reads[role] as List).size()]
        }
        if( counts.values().toSet().size() != 1 ) {
            final String countSummary = roles.collect { role -> "${role}=${counts[role]}" }.join(', ')
            throw new IllegalArgumentException(
                "samples.${sampleId}.${modality}.reads has conflicting technical read-set counts: ${countSummary}"
            )
        }

        final int readSetCount = counts[roles[0]]
        (0..<readSetCount).each { index ->
            final List<Map<String, String>> seenPaths = []
            roles.each { role ->
                if( !(implicitDualI2 && role == 'i2') ) {
                    final String path = (reads[role] as List<String>)[index]
                    final Map<String, String> conflict = seenPaths.find { seen ->
                        java.nio.file.Files.isSameFile(
                            new File(seen.path).toPath(),
                            new File(path).toPath()
                        )
                    }
                    if( conflict ) {
                        throw new IllegalArgumentException(
                            "samples.${sampleId}.${modality}.reads reuses one physical FASTQ in technical " +
                            "read set ${index + 1} for roles ${conflict.role} and ${role}: ${path}"
                        )
                    }
                    seenPaths << [role: role, path: path]
                }
            }
        }
    }

    private static Map resolveRuntime(final Map parsed, final File baseDir, final Map options) {
        final Map runtime = asMap(parsed.runtime, 'runtime')
        final String explicitTmpdir = runtime.tmpdir?.toString()?.trim()
        final String defaultTmpdir = options.outdir?.toString()?.trim() ?: 'results'

        return [
            env_prefix: resolvePath(
                baseDir,
                requireString(runtime.env_prefix, 'runtime.env_prefix')
            ),
            tmpdir: resolvePath(
                baseDir,
                explicitTmpdir ?: defaultTmpdir
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
