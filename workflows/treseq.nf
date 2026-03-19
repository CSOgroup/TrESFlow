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
 */

include { INITIAL_RNA_TAGGING } from '../subworkflows/local/initial_rna_tagging'

workflow TRESEQ {
    main:
    if( !params.samplesheet ) {
        error "Missing required parameter: --samplesheet"
    }

    final String species = (params.rna_align_species ?: '').toString().trim().toLowerCase()
    final String refBaseDir = (params.rna_ref_base_dir ?: '').toString().trim()
    final int maxCpus = params.max_cpus as int

    if( !species ) {
        error "Missing required parameter for RNA alignment: --rna_align_species human|mouse"
    }
    if( !(species in ['human', 'mouse']) ) {
        error "Invalid --rna_align_species '${species}'. Supported values: human, mouse"
    }
    if( !refBaseDir ) {
        error "Missing required parameter for RNA alignment: --rna_ref_base_dir"
    }
    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }

    final List<String> requiredRefPaths = species == 'human'
        ? ['GRCh38_TrES/star', 'hg38.chrom.sizes']
        : ['GRCm39_TrES/star', 'mm39.chrom.sizes']

    requiredRefPaths.each { relPath ->
        final File resolved = new File(refBaseDir, relPath)
        if( !resolved.exists() ) {
            error "RNA alignment reference path not found for species '${species}': ${resolved}"
        }
    }

    def sampleRows = SamplesheetParser.parse(params.samplesheet as String)

    Channel
        .fromList(sampleRows)
        .map { row ->
            tuple(
                row.id,
                row,
                file(row.i1),
                file(row.r1),
                file(row.r2),
                file(row.cell_whitelist),
                file(row.rna_sb_group_map)
            )
        }
        .set { ch_rna_samples }

    INITIAL_RNA_TAGGING(ch_rna_samples)

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
}
