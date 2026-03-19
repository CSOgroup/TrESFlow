import groovy.yaml.YamlSlurper

class SamplesheetParser {

    static List<Map> parse(final String samplesheetPath) {
        if( !samplesheetPath ) {
            throw new IllegalArgumentException("Missing required parameter: --samplesheet")
        }

        final File sheetFile = new File(samplesheetPath)
        if( !sheetFile.exists() ) {
            throw new IllegalArgumentException("Samplesheet not found: ${samplesheetPath}")
        }

        final def parsed = new YamlSlurper().parse(sheetFile)
        if( !(parsed instanceof Map) || !(parsed.samples instanceof List) || parsed.samples.isEmpty() ) {
            throw new IllegalArgumentException(
                "Samplesheet must contain a non-empty top-level 'samples:' list: ${samplesheetPath}"
            )
        }

        final File baseDir = sheetFile.parentFile ?: new File('.')

        final List<Map> samples = []

        for( int idx = 0; idx < parsed.samples.size(); idx++ ) {
            final def row = parsed.samples[idx]
            if( !(row instanceof Map) ) {
                throw new IllegalArgumentException("samples[${idx}] must be a mapping")
            }

            final String sampleId = requireString(row.id, "samples[${idx}].id")
            final String modality = requireString(row.modality, "samples[${idx}].modality").toLowerCase()

            final Map reads = asMap(row.reads, "samples[${idx}].reads")
            final Map barcodes = asMap(row.barcodes, "samples[${idx}].barcodes")
            final Map sampleBarcode = asMap(barcodes.sample, "samples[${idx}].barcodes.sample")
            final Map umiBarcode = asMap(barcodes.umi, "samples[${idx}].barcodes.umi")

            final Map sample = [
                id                         : sampleId,
                modality                   : modality,
                r1                         : resolveExistingPath(baseDir, requireString(reads.r1, "samples[${idx}].reads.r1")),
                r2                         : resolveExistingPath(baseDir, requireString(reads.r2, "samples[${idx}].reads.r2")),
                sample_whitelist           : resolveExistingPath(
                    baseDir,
                    requireString(sampleBarcode.whitelist, "samples[${idx}].barcodes.sample.whitelist")
                ),
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
                umi_tag                    : optionalString(umiBarcode.tag, 'UM')
            ]

            if( sample.modality != 'rna' ) {
                throw new IllegalArgumentException(
                    "This first implementation slice only supports modality: rna (sample '${sampleId}' has '${sample.modality}')"
                )
            }

            samples << sample
        }

        return samples
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
