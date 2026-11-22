import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Path

class PreprocessedManifest {
    static final String NAME = 'restart_manifest.json'

    private static void protect(Path source, List<Path> destinations, String label) {
        Path realSource = source.toRealPath()
        if( destinations.any { source.startsWith(it) || realSource.startsWith(it) } ) {
            throw new IllegalArgumentException("${label}: use a fresh --outdir and -work-dir outside the restart inputs: ${source}")
        }
    }

    private static Path resolveInput(Path base, Object value, String label, List<Path> destinations) {
        if( !(value instanceof String) || !value.trim() || Path.of(value).isAbsolute() ) {
            throw new IllegalArgumentException("${label}: expected a path relative to the restart manifest")
        }
        Path source = base.resolve(value).normalize()
        if( !Files.isRegularFile(source) || !Files.isReadable(source) ) {
            throw new IllegalArgumentException("${label}: file is missing or unreadable: ${source}")
        }
        protect(source, destinations, label)
        // Keep the listed basename and compression suffix even when it is a symlink.
        return source
    }

    static Map read(String input, String output, String workDir) {
        Path requested = Path.of(input)
        List<Path> destinations = [Path.of(output), new File(workDir).canonicalFile.toPath()]
        Path manifest = Files.isDirectory(requested) ? requested.resolve(NAME) : requested
        if( !Files.isRegularFile(manifest) ) {
            throw new IllegalArgumentException("Restart manifest not found: ${manifest}. Use a completed preprocessing output or a version 1 manifest.")
        }
        protect(manifest, destinations, 'Restart manifest')
        def doc
        try {
            doc = new JsonSlurper().parse(manifest.toFile())
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Cannot read restart manifest ${manifest}: ${e.message}", e)
        }
        if( !(doc instanceof Map) || doc.version != 1 || !(doc.splits instanceof List) || !doc.splits ) {
            throw new IllegalArgumentException("Restart manifest ${manifest} requires version 1 and a nonempty splits list.")
        }
        if( doc.containsKey('filter_dual_tag_artifacts') && !(doc.filter_dual_tag_artifacts instanceof Boolean) ) {
            throw new IllegalArgumentException("Restart manifest ${manifest}: filter_dual_tag_artifacts must be a boolean.")
        }

        Set readFiles = new HashSet()
        Path base = manifest.parent
        List splits = doc.splits.collect { entry ->
            if( !(entry instanceof Map) ) throw new IllegalArgumentException('Each restart split must be an object.')
            String key = "${entry.modality}:${entry.split_id}"
            String label = "Restart split '${key}'"
            if( entry.containsKey('pairs') && (!(entry.pairs instanceof Integer || entry.pairs instanceof Long) || entry.pairs < 0) ) {
                throw new IllegalArgumentException("${label}: pairs must be a nonnegative integer.")
            }
            def reads = ['r1', 'r2'].collect { role ->
                Path source = resolveInput(base, entry[role], "${label} ${role}", destinations)
                if( !(source.fileName.toString() ==~ /.+\.(?:fastq|fq)(?:\.gz)?/) ) {
                    throw new IllegalArgumentException("${label} ${role}: expected .fastq or .fq, optionally gzip compressed: ${source}")
                }
                def physicalFile = Files.readAttributes(source, java.nio.file.attribute.BasicFileAttributes).fileKey() ?: source.toRealPath()
                if( !readFiles.add(physicalFile) ) {
                    throw new IllegalArgumentException("${label} ${role}: FASTQ is reused in this manifest: ${source}")
                }
                source
            }
            Path rg = entry.rg == null ? null : resolveInput(base, entry.rg, "${label} rg", destinations)
            [entry: entry, r1: reads[0], r2: reads[1], rg: rg]
        }
        Map qc = ['report_inputs', 'multiqc_inputs'].collectEntries { key ->
            def values = doc.containsKey(key) ? doc[key] : []
            if( !(values instanceof List) || values.toSet().size() != values.size() ) {
                throw new IllegalArgumentException("Restart manifest ${manifest}: ${key} must be a list of distinct paths.")
            }
            [(key): values.collect { resolveInput(base, it, "Restart ${key}", destinations) }]
        }
        return [splits: splits, filter_dual_tag_artifacts: doc.filter_dual_tag_artifacts] + qc
    }

    static void bindSamples(Map restart, List<Map> rows) {
        Map expected = rows.collectMany { row ->
            row.split_targets.collect { id, target -> [row.modality + ':' + id, [row: row, target: target]] }
        }.collectEntries()
        restart.splits.each { split ->
            Map entry = split.entry
            String key = "${entry.modality}:${entry.split_id}"
            String label = "Restart split '${key}'"
            Map spec = expected.remove(key)
            if( !spec ) {
                throw new IllegalArgumentException("${label}: duplicate or absent from the samplesheet split definitions.")
            }
            Map row = spec.row
            Map target = spec.target
            boolean dna = row.modality == 'dna'
            Map identity = [
                sample: row.id, library: row.library_name, group: target.group,
                format: 'tresflow-split-v1',
            ]
            if( dna ) identity += [mark: target.mark, dna_tagmentation: row.dna_tagmentation]
            identity.each { field, value ->
                if( entry[field] != value ) {
                    throw new IllegalArgumentException("${label}: expected ${field} '${value}', received '${entry[field]}'. Supply exact split FASTQs and matching samplesheet metadata.")
                }
            }
            if( !dna && (entry.mark != null || entry.rg != null || entry.dna_tagmentation != null) ) {
                throw new IllegalArgumentException("${label}: mark, dna_tagmentation and rg are DNA-only fields.")
            }
            if( entry.containsKey('sb_mapping_sha256') && entry.sb_mapping_sha256 != row.sb_mapping_sha256 ) {
                throw new IllegalArgumentException("${label}: samplesheet SB mapping differs from preprocessing.")
            }
            split.row = row
        }
        if( expected ) {
            throw new IllegalArgumentException("Restart manifest is missing samplesheet splits: ${expected.keySet().join(', ')}. Match the samplesheet to the splits you intend to restart.")
        }
    }
}
