/*
 * Subworkflow: OPTIONAL_SC_PROCESS
 * Optional boundary:
 *   - Runs only when shared downstream staging or the optional sc_process path is enabled.
 *
 * Inputs:
 *   - RNA grouped STARsolo directories
 *   - RNA grouped filtered BAMs
 *   - DNA grouped NoDup BAMs
 *   - shared metadata: label, species, genome, mo_map, sb_group_map
 *   - boolean controlling stage-only vs stage+run behavior
 * Outputs:
 *   - one flat shared stage directory
 *   - one optional sc_process working directory when --run_sc_process is enabled
 */

import WorkflowSupport

include { STAGE_SC_PROCESS_INPUTS } from '../../modules/local/stage_sc_process_inputs/main'
include { RUN_SC_PROCESS } from '../../modules/local/run_sc_process/main'

workflow OPTIONAL_SC_PROCESS {
    take:
    ch_rna_solo_dirs
    ch_rna_filtered_bams
    ch_dna_nodup_bams
    ch_stage_meta
    runScProcess

    main:
    ch_rna_solo_bundle = ch_rna_solo_dirs
        .map { splitName, meta, soloDir -> soloDir }
        .collect()
        .map { soloDirs -> [rnaSoloDirs: WorkflowSupport.asPathList(soloDirs)] }

    ch_rna_filtered_bundle = ch_rna_filtered_bams
        .map { splitName, meta, filteredBam -> filteredBam }
        .collect()
        .map { filteredBams -> [rnaFilteredBams: WorkflowSupport.asPathList(filteredBams)] }

    ch_dna_nodup_bundle = ch_dna_nodup_bams
        .map { splitName, meta, noDupBam -> noDupBam }
        .collect()
        .map { noDupBams -> [dnaNoDupBams: WorkflowSupport.asPathList(noDupBams)] }

    ch_stage_input = ch_stage_meta
        .combine(ch_rna_solo_bundle)
        .combine(ch_rna_filtered_bundle)
        .combine(ch_dna_nodup_bundle)
        .map { values ->
            tuple(
                values[0],
                values[1],
                values[2],
                values[3],
                values[4],
                values[5].rnaSoloDirs,
                values[6].rnaFilteredBams,
                values[7].dnaNoDupBams
            )
        }

    STAGE_SC_PROCESS_INPUTS(ch_stage_input)

    def ch_run_dir = Channel.empty()
    if( runScProcess ) {
        ch_sc_process_input = ch_stage_meta
            .combine(STAGE_SC_PROCESS_INPUTS.out.stage_dir)
            .map { values ->
                tuple(
                    values[0],
                    values[1],
                    values[2],
                    values[5]
                )
            }

        RUN_SC_PROCESS(ch_sc_process_input)
        ch_run_dir = RUN_SC_PROCESS.out.run_dir
    }

    emit:
    stage_dir = STAGE_SC_PROCESS_INPUTS.out.stage_dir
    run_dir = ch_run_dir
}
