/*
 * Workflow: TRESEQ
 * Current slice:
 *   1. Parse a single YAML samplesheet.
 *   2. Run the upstream RNA sample-barcode tagging step (Tag.codon) via a thin wrapper.
 *   3. Run the upstream RNA UMI tagging step (Tag_UMI.codon) via a thin wrapper.
 *   4. Run the upstream RNA cell-barcode tagging step (Tag_Lig3.codon) via a thin wrapper.
 *   5. Run the upstream RNA trim_galore step via a thin wrapper.
 *   6. Run the upstream RNA Split_ReadsV2 step in rna mode via a thin wrapper.
 *   7. Run the upstream RNA FqToSAM step via a thin wrapper.
 *   8. Run the upstream RNA AlignRNA.sh step via a thin wrapper.
 *   9. Run the upstream DNA sample-barcode, modality-barcode, and cell-barcode tagging
 *      steps plus DNA trim_galore, Split_ReadsV2 dna mode, AlignDNA.sh,
 *      and GATK MarkDuplicates through the duplicate-marked DNA BAM boundary.
 */

import PipelineSupport

include { INITIAL_RNA_TAGGING } from '../subworkflows/local/initial_rna_tagging'
include { INITIAL_DNA_TAGGING } from '../subworkflows/local/initial_dna_tagging'

def toRnaSampleInput(final Map row) {
    tuple(
        row.id,
        row,
        file(row.i1),
        file(row.r1),
        file(row.r2),
        file(row.cell_whitelist),
        file(row.sb_group_map)
    )
}

def toDnaSampleInput(final Map row) {
    tuple(
        row.id,
        row,
        file(row.i1),
        file(row.i2),
        file(row.r1),
        file(row.r2),
        file(row.modality_whitelist),
        file(row.cell_whitelist),
        file(row.mo_map),
        file(row.sb_group_map)
    )
}

workflow TRESEQ {
    main:
    if( !params.samplesheet ) {
        error "Missing required parameter: --samplesheet"
    }

    def sampleRows = SamplesheetParser.parse(params.samplesheet as String)
    final List<Map> rnaRows = sampleRows.findAll { row -> row.modality == 'rna' }
    final List<Map> dnaRows = sampleRows.findAll { row -> row.modality == 'dna' }
    final int maxCpus = params.max_cpus as int

    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }

    if( rnaRows ) {
        PipelineSupport.validateRnaAlignment(
            params.rna_ref_base_dir as String,
            params.rna_align_species as String
        )
    }

    if( dnaRows ) {
        try {
            PipelineSupport.validateDnaAlignment(
                params.dna_bwa_reference as String,
                params.dna_blacklist_bed as String,
                params.dna_effective_genome_size as String
            )
        }
        catch( IllegalArgumentException e ) {
            error e.message
        }
    }

    Channel
        .fromList(rnaRows)
        .map { row -> toRnaSampleInput(row) }
        .set { ch_rna_samples }

    Channel
        .fromList(dnaRows)
        .map { row -> toDnaSampleInput(row) }
        .set { ch_dna_samples }

    INITIAL_RNA_TAGGING(ch_rna_samples)
    INITIAL_DNA_TAGGING(ch_dna_samples)

    emit:
    tagged_fastqs     = INITIAL_RNA_TAGGING.out.tagged_fastqs
    trimmed_fastqs    = INITIAL_RNA_TAGGING.out.trimmed_fastqs
    split_fastqs      = INITIAL_RNA_TAGGING.out.split_fastqs
    rg_headers        = INITIAL_RNA_TAGGING.out.rg_headers
    usam_files        = INITIAL_RNA_TAGGING.out.usam_files
    aligned_solo_dirs = INITIAL_RNA_TAGGING.out.aligned_solo_dirs
    aligned_filtered_bams = INITIAL_RNA_TAGGING.out.aligned_filtered_bams
    aligned_stranded_bigwigs = INITIAL_RNA_TAGGING.out.aligned_stranded_bigwigs
    aligned_unstranded_bigwigs = INITIAL_RNA_TAGGING.out.aligned_unstranded_bigwigs
    barcode_reports   = INITIAL_RNA_TAGGING.out.barcode_reports
    dna_tagged_fastqs = INITIAL_DNA_TAGGING.out.tagged_fastqs
    dna_trimmed_fastqs = INITIAL_DNA_TAGGING.out.trimmed_fastqs
    dna_split_fastqs = INITIAL_DNA_TAGGING.out.split_fastqs
    dna_rg_headers = INITIAL_DNA_TAGGING.out.rg_headers
    dna_aligned_bams = INITIAL_DNA_TAGGING.out.aligned_bams
    dna_aligned_bais = INITIAL_DNA_TAGGING.out.aligned_bais
    dna_alignment_barcode_counts = INITIAL_DNA_TAGGING.out.alignment_barcode_counts
    dna_markeddup_bams = INITIAL_DNA_TAGGING.out.markeddup_bams
    dna_markeddup_bais = INITIAL_DNA_TAGGING.out.markeddup_bais
    dna_duplicate_metrics = INITIAL_DNA_TAGGING.out.duplicate_metrics
    dna_barcode_reports = INITIAL_DNA_TAGGING.out.barcode_reports
}
