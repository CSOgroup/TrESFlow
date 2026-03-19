/*
 * Workflow: TRESEQ
 * Current slice:
 *   1. Parse a single YAML samplesheet.
 *   2. Run the upstream RNA sample-barcode tagging step (Tag.codon) via a thin wrapper.
 *   3. Run the upstream RNA UMI tagging step (Tag_UMI.codon) via a thin wrapper.
 *   4. Run the upstream RNA cell-barcode tagging step (Tag_Lig3.codon) via a thin wrapper.
 *   5. Run the upstream RNA trim_galore step via a thin wrapper.
 */

include { INITIAL_RNA_TAGGING } from '../subworkflows/local/initial_rna_tagging'

workflow TRESEQ {
    main:
    if( !params.samplesheet ) {
        error "Missing required parameter: --samplesheet"
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
                file(row.sample_whitelist),
                file(row.cell_whitelist)
            )
        }
        .set { ch_rna_samples }

    INITIAL_RNA_TAGGING(ch_rna_samples)

    emit:
    tagged_fastqs     = INITIAL_RNA_TAGGING.out.tagged_fastqs
    trimmed_fastqs    = INITIAL_RNA_TAGGING.out.trimmed_fastqs
    barcode_reports   = INITIAL_RNA_TAGGING.out.barcode_reports
}
