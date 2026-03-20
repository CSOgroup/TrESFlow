# Implemented Pipeline Architecture

Repo-maintained architecture view for the current implementation.

- Core workflow: RNA through `ALIGN_RNA`, DNA through `BAM_COVERAGE_DNA`
- Optional downstream: shared staging plus one optional `sc_process.py` path

```mermaid
flowchart TD
    SS[Single YAML samplesheet]
    GROUPS[groups with sb_barcodes]
    DNAMARKS[DNA mark_barcodes]
    DERIVE[Derive internal contract files\nsb_group_map.tsv\ndna_mo_map.tsv\nDNA modality whitelists]
    LIG[ligation_barcode_whitelist]
    RNAREF[rna_ref_base_dir]
    DNAREF[dna_bwa_reference]
    DNABL[dna_blacklist_bed]
    DNAEFF[dna_effective_genome_size]

    SS --> GROUPS
    SS --> RNA0
    SS --> DNA0
    SS --> DNAMARKS
    GROUPS --> DERIVE
    DNAMARKS --> DERIVE
    DERIVE --> RNA0
    DERIVE --> RNA4
    DERIVE --> DNA0
    DERIVE --> DNA1
    DERIVE --> DNA4
    LIG --> RNA2
    LIG --> DNA2
    RNAREF --> RNA6
    DNAREF --> DNA5
    DNABL --> DNA5
    DNAEFF --> DNA5
    DNAEFF --> DNA8

    subgraph RNA_Core[RNA Core Branch]
        RNA0[TAG_RNA_SAMPLE_BARCODE]
        RNA1[TAG_RNA_UMI]
        RNA2[TAG_RNA_CELL_BARCODE]
        RNA3[TRIM_RNA_FASTQS]
        RNA4[SPLIT_RNA_READS]
        RNA5[FQ_TO_SAM]
        RNA6[ALIGN_RNA]
        RNA0 --> RNA1 --> RNA2 --> RNA3 --> RNA4 --> RNA5 --> RNA6
    end

    subgraph DNA_Core[DNA Core Branch]
        DNA0[TAG_DNA_SAMPLE_BARCODE]
        DNA1[TAG_DNA_MODALITY_BARCODE]
        DNA2[TAG_DNA_CELL_BARCODE]
        DNA3[TRIM_DNA_FASTQS]
        DNA4[SPLIT_DNA_READS]
        DNA5[ALIGN_DNA]
        DNA6[MARK_DUPLICATES_DNA]
        DNA7[SPLIT_DUPLICATES_DNA]
        DNA8[BAM_COVERAGE_DNA]
        DNA0 --> DNA1 --> DNA2 --> DNA3 --> DNA4 --> DNA5 --> DNA6 --> DNA7 --> DNA8
    end

    RNA6 --> STAGE
    DNA7 --> STAGE

    subgraph Optional_Downstream[Optional Downstream]
        STAGE[STAGE_SC_PROCESS_INPUTS]
        RUNSC[RUN_SC_PROCESS]
        STAGE --> RUNSC
    end
```

Notes:

- RNA and DNA stay independent and parallel until the optional shared downstream boundary.
- The core workflow does not require `sc_process.py`.
- `sb_group_map.tsv`, `dna_mo_map.tsv`, and DNA modality whitelist files are derived internally from the hierarchical YAML.
- The core runtime scripts are repo-owned under `scripts/core_runtime/`.
- The optional downstream `sc_process.py` path remains separate by design.
