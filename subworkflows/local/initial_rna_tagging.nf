/*
 * Subworkflow: INITIAL_RNA_TAGGING
 * Inputs:
 *   - sample metadata parsed from params.samplesheet
 *   - raw RNA I1 / R1 / R2 FASTQs
 *   - RNA cell-barcode whitelist
 *   - RNA SB-group map TSV, used both to derive the effective sample-barcode whitelist and later by Split_ReadsV2 rna mode
 * Outputs:
 *   - RNA FASTQs tagged with SB, UM, then CB comments
 *   - trim_galore paired-end FASTQs from the CB-tagged reads
 *   - Split_ReadsV2 per-group RNA FASTQs and SAM RG headers
 *   - FqToSAM unmapped SAM files from each split RNA FASTQ pair
 *   - barcode count/stat files from all wrapped upstream steps
 */

include { TAG_RNA_SAMPLE_BARCODE } from '../../modules/local/tag_rna_sb/main'
include { TAG_RNA_UMI }            from '../../modules/local/tag_rna_umi/main'
include { TAG_RNA_CELL_BARCODE }   from '../../modules/local/tag_rna_cell_barcode/main'
include { TRIM_RNA_FASTQS }        from '../../modules/local/trim_rna_fastqs/main'
include { SPLIT_RNA_READS }        from '../../modules/local/split_rna_reads/main'
include { FQ_TO_SAM }              from '../../modules/local/fq_to_sam/main'

def asPathList(obj) {
    if( obj instanceof List ) {
        return obj
    }
    return [obj]
}

workflow INITIAL_RNA_TAGGING {
    take:
    ch_rna_samples

    main:
    ch_sb_input = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, r1, r2, sbGroupMap)
    }

    TAG_RNA_SAMPLE_BARCODE(ch_sb_input)

    ch_raw_r2 = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, r2)
    }

    ch_umi_input = ch_raw_r2
        .join(TAG_RNA_SAMPLE_BARCODE.out.tagged)
        .map { sampleId, metaFromInput, rawR2, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, rawR2, taggedR1, taggedR2)
        }

    TAG_RNA_UMI(ch_umi_input)

    ch_cb_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, i1, cellWhitelist)
    }

    ch_cb_input = ch_cb_meta
        .join(TAG_RNA_UMI.out.tagged)
        .map { sampleId, metaFromInput, i1, cellWhitelist, metaFromTag, taggedR1, taggedR2 ->
            tuple(sampleId, metaFromInput, i1, taggedR1, taggedR2, cellWhitelist)
        }

    TAG_RNA_CELL_BARCODE(ch_cb_input)
    TRIM_RNA_FASTQS(TAG_RNA_CELL_BARCODE.out.tagged)

    ch_split_meta = ch_rna_samples.map { sampleId, meta, i1, r1, r2, cellWhitelist, sbGroupMap ->
        tuple(sampleId, meta, sbGroupMap)
    }

    ch_split_input = ch_split_meta
        .join(TRIM_RNA_FASTQS.out.trimmed)
        .map { sampleId, metaFromInput, sbGroupMap, metaFromTrim, trimmedR1, trimmedR2 ->
            tuple(sampleId, metaFromInput, trimmedR1, trimmedR2, sbGroupMap)
        }

    SPLIT_RNA_READS(ch_split_input)

    ch_fq_to_sam_input = SPLIT_RNA_READS.out.split_fastqs
        .flatMap { sampleId, meta, splitR1s, splitR2s ->
            def r1ByGroup = asPathList(splitR1s).collectEntries { path ->
                def name = path.getName()
                def group = name.replaceFirst("^${sampleId}_", '').replaceFirst('_R1\\.fq\\.gz$', '')
                [(group): path]
            }
            def r2ByGroup = asPathList(splitR2s).collectEntries { path ->
                def name = path.getName()
                def group = name.replaceFirst("^${sampleId}_", '').replaceFirst('_R2\\.fq\\.gz$', '')
                [(group): path]
            }

            def groups = (r1ByGroup.keySet() + r2ByGroup.keySet()).unique().sort()
            groups.collect { group ->
                if( !r1ByGroup.containsKey(group) || !r2ByGroup.containsKey(group) ) {
                    throw new IllegalStateException("Missing split FASTQ mate for sample '${sampleId}' group '${group}'")
                }
                tuple("${sampleId}_${group}", meta, r1ByGroup[group], r2ByGroup[group])
            }
        }

    FQ_TO_SAM(ch_fq_to_sam_input)

    ch_barcode_reports = TAG_RNA_SAMPLE_BARCODE.out.metrics
        .mix(TAG_RNA_UMI.out.metrics)
        .mix(TAG_RNA_CELL_BARCODE.out.metrics)

    emit:
    tagged_fastqs    = TAG_RNA_CELL_BARCODE.out.tagged
    trimmed_fastqs   = TRIM_RNA_FASTQS.out.trimmed
    split_fastqs     = SPLIT_RNA_READS.out.split_fastqs
    rg_headers       = SPLIT_RNA_READS.out.rg_headers
    usam_files       = FQ_TO_SAM.out.usam
    barcode_reports  = ch_barcode_reports
}
