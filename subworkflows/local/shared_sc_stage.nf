/*
 * Subworkflow: SHARED_SC_STAGE
 * Optional boundary:
 *   - Runs only when the user explicitly enables shared downstream staging.
 *
 * Inputs:
 *   - RNA grouped STARsolo directories
 *   - RNA grouped filtered BAMs
 *   - DNA grouped NoDup BAMs
 *   - shared stage metadata: label, species, genome, mo_map, sb_group_map
 * Outputs:
 *   - one flat shared stage directory ready for a future single sc_process.py invocation
 */

import PipelineSupport

include { STAGE_SC_PROCESS_INPUTS } from '../../modules/local/stage_sc_process_inputs/main'

workflow SHARED_SC_STAGE {
    take:
    ch_rna_solo_dirs
    ch_rna_filtered_bams
    ch_dna_nodup_bams
    ch_stage_meta

    main:
    ch_rna_solo_bundle = ch_rna_solo_dirs
        .map { splitName, meta, soloDir -> soloDir }
        .collect()
        .map { soloDirs -> [rnaSoloDirs: PipelineSupport.asPathList(soloDirs)] }

    ch_rna_filtered_bundle = ch_rna_filtered_bams
        .map { splitName, meta, filteredBam -> filteredBam }
        .collect()
        .map { filteredBams -> [rnaFilteredBams: PipelineSupport.asPathList(filteredBams)] }

    ch_dna_nodup_bundle = ch_dna_nodup_bams
        .map { splitName, meta, noDupBam -> noDupBam }
        .collect()
        .map { noDupBams -> [dnaNoDupBams: PipelineSupport.asPathList(noDupBams)] }

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

    emit:
    stage_dir = STAGE_SC_PROCESS_INPUTS.out.stage_dir
}
