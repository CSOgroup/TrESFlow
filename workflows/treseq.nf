/*
 * Workflow: TRESEQ
 * Current slice:
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
 *      GATK MarkDuplicates, duplicate filtering to NoDup BAMs, and bamCoverage
 *      through the DNA-only post-alignment prep boundary before shared downstream work.
 *  10. Optionally stage one flat shared downstream workdir.
 *  11. Optionally run one and only one sc_process.py call after that shared stage exists.
 */

import WorkflowSupport
import RuntimeSupport

include { RNA_CORE } from '../subworkflows/local/rna_core'
include { DNA_CORE } from '../subworkflows/local/dna_core'
include { OPTIONAL_SC_PROCESS } from '../subworkflows/local/optional_sc_process'

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

def parseBooleanParam(final Object rawValue, final String paramName) {
    if( rawValue instanceof Boolean ) {
        return (Boolean) rawValue
    }

    final String normalized = rawValue?.toString()?.trim()?.toLowerCase()
    if( !normalized ) {
        return false
    }
    if( normalized in ['true', '1', 'yes', 'y'] ) {
        return true
    }
    if( normalized in ['false', '0', 'no', 'n'] ) {
        return false
    }

    error "Invalid --${paramName} '${rawValue}'. Supported boolean values: true, false"
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
        ]
    )
    final List<Map> rnaRows = sampleRows.findAll { row -> row.modality == 'rna' }
    final List<Map> dnaRows = sampleRows.findAll { row -> row.modality == 'dna' }
    final int maxCpus = params.max_cpus as int
    final boolean stageScProcessInputs = parseBooleanParam(params.stage_sc_process_inputs, 'stage_sc_process_inputs')
    final boolean runScProcess = parseBooleanParam(params.run_sc_process, 'run_sc_process')
    final boolean enableSharedStage = stageScProcessInputs || runScProcess

    if( maxCpus < 1 ) {
        error "Invalid --max_cpus '${maxCpus}'. Value must be >= 1"
    }

    if( enableSharedStage && !(rnaRows && dnaRows) ) {
        final String flagName = runScProcess ? '--run_sc_process' : '--stage_sc_process_inputs'
        error "${flagName} requires both RNA and DNA samples in the same YAML samplesheet"
    }

    if( rnaRows ) {
        WorkflowSupport.validateRnaAlignment(
            params.rna_ref_base_dir as String,
            params.rna_align_species as String
        )
    }

    if( dnaRows ) {
        try {
            WorkflowSupport.validateDnaAlignment(
                params.dna_bwa_reference as String,
                params.dna_blacklist_bed as String,
                params.dna_effective_genome_size as String
            )
        }
        catch( IllegalArgumentException e ) {
            error e.message
        }
    }

    if( runScProcess ) {
        try {
            RuntimeSupport.validateConfiguredDirectory(
                'runtime SnapATAC cache',
                params.runtime_snap_data_dir as String
            )
        }
        catch( IllegalStateException e ) {
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

    def ch_shared_stage_dir = Channel.empty()
    def ch_sc_process_run_dir = Channel.empty()
    if( enableSharedStage ) {
        final String sharedSpecies = (params.rna_align_species ?: '').toString().trim().toLowerCase()
        final String sharedGenome = WorkflowSupport.sharedGenomeFromSpecies(sharedSpecies)
        final String sharedStageLabel = (sampleRows[0].library_name ?: 'shared_stage').toString()
        final String sharedSbGroupMap = (dnaRows ? dnaRows[0].sb_group_map : rnaRows[0].sb_group_map).toString()
        final String sharedMoMap = dnaRows[0].mo_map.toString()

        Channel
            .value(
                tuple(
                    sharedStageLabel,
                    sharedSpecies,
                    sharedGenome,
                    file(sharedMoMap),
                    file(sharedSbGroupMap)
                )
            )
            .set { ch_shared_stage_meta }

        OPTIONAL_SC_PROCESS(
            RNA_CORE.out.aligned_solo_dirs,
            RNA_CORE.out.aligned_filtered_bams,
            DNA_CORE.out.nodup_bams,
            ch_shared_stage_meta,
            runScProcess
        )
        ch_shared_stage_dir = OPTIONAL_SC_PROCESS.out.stage_dir
        ch_sc_process_run_dir = OPTIONAL_SC_PROCESS.out.run_dir
    }

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
    shared_sc_stage = ch_shared_stage_dir
    sc_process_run = ch_sc_process_run_dir
}
