/*
 * Workflow: TRESEQ
 * Core workflow:
 *   1. Parse a single hierarchical YAML samplesheet.
 *   2. Run the upstream RNA sample-barcode tagging step (Tag.codon) via a thin wrapper.
 *   3. Run the upstream RNA UMI tagging step (Tag_UMI.codon) via a thin wrapper.
 *   4. Run the upstream RNA cell-barcode tagging step (Tag_Lig3.codon) via a thin wrapper.
 *   5. Run the upstream RNA trim_galore step via a thin wrapper.
 *   6. Run the upstream RNA Split_ReadsV2 step in rna mode via a thin wrapper.
 *   7. Run the upstream RNA FqToSAM step via a thin wrapper.
 *   8. Run the upstream RNA AlignRNA.sh step via a thin wrapper.
 *   9. Run the upstream DNA sample-barcode, modality-barcode, and cell-barcode tagging
 *      steps plus DNA trim_galore, Split_ReadsV2 dna mode, AlignDNA.sh,
 *      GATK MarkDuplicates, duplicate filtering to NoDup BAMs, and bamCoverage.
 */

import WorkflowSupport

include { RNA_CORE } from '../subworkflows/local/rna_core'
include { DNA_CORE } from '../subworkflows/local/dna_core'

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

    def sampleRows = SamplesheetParser.parse(
        params.samplesheet as String,
        [
            outdir                     : params.outdir,
            ligation_barcode_whitelist : params.ligation_barcode_whitelist,
            barcode_defaults           : params.barcode_defaults,
            rna_ref_base_dir           : params.rna_ref_base_dir,
            rna_align_species          : params.rna_align_species,
            dna_bwa_reference          : params.dna_bwa_reference,
            dna_blacklist_bed          : params.dna_blacklist_bed,
            dna_effective_genome_size  : params.dna_effective_genome_size,
        ]
    )
    final List<Map> rnaRows = sampleRows.findAll { row -> row.modality == 'rna' }
    final List<Map> dnaRows = sampleRows.findAll { row -> row.modality == 'dna' }
    final int maxCpus = params.max_cpus as int

    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }

    if( rnaRows ) {
        WorkflowSupport.validateRnaAlignment(
            rnaRows[0].rna_ref_base_dir as String,
            rnaRows[0].rna_align_species as String
        )
    }

    if( dnaRows ) {
        try {
            WorkflowSupport.validateDnaAlignment(
                dnaRows[0].dna_bwa_reference as String,
                dnaRows[0].dna_blacklist_bed as String,
                dnaRows[0].dna_effective_genome_size as String
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

    RNA_CORE(ch_rna_samples)
    DNA_CORE(ch_dna_samples)

    emit:
    tagged_fastqs     = RNA_CORE.out.tagged_fastqs
    trimmed_fastqs    = RNA_CORE.out.trimmed_fastqs
    split_fastqs      = RNA_CORE.out.split_fastqs
    rg_headers        = RNA_CORE.out.rg_headers
    usam_files        = RNA_CORE.out.usam_files
    aligned_solo_dirs = RNA_CORE.out.aligned_solo_dirs
    aligned_filtered_bams = RNA_CORE.out.aligned_filtered_bams
    aligned_stranded_bigwigs = RNA_CORE.out.aligned_stranded_bigwigs
    aligned_unstranded_bigwigs = RNA_CORE.out.aligned_unstranded_bigwigs
    barcode_reports   = RNA_CORE.out.barcode_reports
    dna_tagged_fastqs = DNA_CORE.out.tagged_fastqs
    dna_trimmed_fastqs = DNA_CORE.out.trimmed_fastqs
    dna_split_fastqs = DNA_CORE.out.split_fastqs
    dna_rg_headers = DNA_CORE.out.rg_headers
    dna_aligned_bams = DNA_CORE.out.aligned_bams
    dna_aligned_bais = DNA_CORE.out.aligned_bais
    dna_alignment_barcode_counts = DNA_CORE.out.alignment_barcode_counts
    dna_markeddup_bams = DNA_CORE.out.markeddup_bams
    dna_markeddup_bais = DNA_CORE.out.markeddup_bais
    dna_duplicate_metrics = DNA_CORE.out.duplicate_metrics
    dna_nodup_bams = DNA_CORE.out.nodup_bams
    dna_nodup_bais = DNA_CORE.out.nodup_bais
    dna_coverage_bigwigs = DNA_CORE.out.coverage_bigwigs
    dna_barcode_reports = DNA_CORE.out.barcode_reports
}
