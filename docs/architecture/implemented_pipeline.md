# Implemented Pipeline Architecture

Core workflow only:

- RNA through the repo-owned STARsolo, filtered-BAM, and coverage stages
- DNA through repo-owned tagging/splitting/alignment, nf-core `gatk4/markduplicates`, repo-owned NoDup extraction, and nf-core `deeptools/bamcoverage`
- nf-core FastQC, combined Samtools sidecar QC, nf-core MultiQC, and a TrESFlow-specific HTML report at the end of the run

```mermaid
flowchart TD
    SS[Hierarchical YAML samplesheet]
    RUN[runtime block]
    REF[references block\nspecies and direct paths]
    GROUPS[groups with sb_barcodes]
    DNAMARKS[DNA mark_barcodes]
    DERIVE[Derive internal files\nbarcode maps and whitelists\ncanonical chromosome contracts]

    SS --> RUN
    SS --> REF
    SS --> GROUPS
    SS --> DNAMARKS
    GROUPS --> DERIVE
    DNAMARKS --> DERIVE

    REF --> RNA2
    REF --> RNA6
    REF --> DNA2
    REF --> DNA5
    REF --> DNA8
    DERIVE --> RNA0
    DERIVE --> RNA4
    DERIVE --> DNA0
    DERIVE --> DNA1
    DERIVE --> DNA4

    subgraph RNA_Core[RNA Core]
        RNA0[TAG_RNA_SAMPLE_BARCODE]
        RNA1[TAG_RNA_UMI]
        RNA2[TAG_RNA_CELL_BARCODE]
        RNA3[TRIM_RNA_FASTQS]
        RNA4[SPLIT_RNA_READS]
        RNA4P[COMPRESS_RNA_SPLIT_FASTQS\noptional published .fastq.gz]
        RNA5[FQ_TO_SAM]
        RNA6[RNA_STARSOLO_ALIGN]
        RNA7[RNA_FILTERED_BAM\ncanonical low-compression published BAM]
        RNA8[RNA_COVERAGE]
        RNA0 --> RNA1 --> RNA2 --> RNA3 --> RNA4 --> RNA5 --> RNA6 --> RNA7 --> RNA8
        RNA4 --> RNA4P
    end

    subgraph DNA_Core[DNA Core]
        DNA0[TAG_DNA_SAMPLE_BARCODE]
        DNA1[TAG_DNA_MODALITY_BARCODE]
        DNA2[TAG_DNA_CELL_BARCODE]
        DNA3[TRIM_DNA_FASTQS]
        DNA4[SPLIT_DNA_READS]
        DNA4P[COMPRESS_DNA_SPLIT_FASTQS\noptional published .fastq.gz]
        DNA5[ALIGN_DNA]
        DNA5C[FILTER_CANONICAL_DNA_ALIGNED_BAM\nQC/output copy]
        DNA6[nf-core GATK4_MARKDUPLICATES\n+ canonical filename normalization]
        DNA7[SPLIT_DUPLICATES_DNA\nNoDup + mapped-read gate]
        DNA8[nf-core DEEPTOOLS_BAMCOVERAGE\ndirect *_NoDup.bw]
        DNA0 --> DNA1 --> DNA2 --> DNA3 --> DNA4 --> DNA5 --> DNA6 --> DNA7 --> DNA8
        DNA5 --> DNA5C
        DNA4 --> DNA4P
    end

    subgraph Reporting[Shared Reporting]
        FASTQC[nf-core FASTQC]
        SAMTOOLS[SAMTOOLS_BAM_QC\none combined task per BAM]
        MULTIQC[nf-core MULTIQC]
        TRESHTML[TRES_REPORT_HTML]
    end

    RNA0 --> FASTQC
    RNA7 --> SAMTOOLS
    RNA6 --> MULTIQC
    DNA0 --> FASTQC
    DNA5C --> SAMTOOLS
    DNA6 --> SAMTOOLS
    DNA7 --> SAMTOOLS
    FASTQC --> MULTIQC
    SAMTOOLS --> MULTIQC
    SAMTOOLS --> TRESHTML
    RNA6 --> TRESHTML
    RNA2 --> TRESHTML
    DNA2 --> TRESHTML
```

Notes:

- One hierarchical samplesheet can describe RNA-only, DNA-only, or combined runs.
- `sb_group_map.tsv`, `dna_mo_map.tsv`, and DNA modality whitelist files are internal artifacts, not user inputs.
- `TAG_DNA_CELL_BARCODE` uses DNA `i1` as the ligation source: single reads use starts `15,53,91`; dual reads use starts `41,79,117`.
- Each Codon split task emits plain FASTQs directly to computation (`FQ_TO_SAM` or `ALIGN_DNA`). With `--publish_split_fastqs`, an independent `pigz -c` branch publishes gzip copies; by default that publication-only task is skipped.
- `RNA_FILTERED_BAM` publishes its low-compression filtered BAM directly under `rna_align/`; RNA coverage and Samtools QC consume the same process output without a copy or recompression task.
- Canonical chromosome contracts are resolved once from the STAR and bwa-mem2 dictionaries. No contigs are renamed, and mixed UCSC/Ensembl conventions fail before task execution.
- GATK duplicate marking consumes the same unfiltered aligned BAM as before. `NORMALIZE_DNA_MARKDUPLICATES` performs the post-marking canonical filter; `SPLIT_DUPLICATES_DNA` therefore removes only records carrying flag `0x400`, indexes the NoDup BAM, and gates coverage from its mapped-read count without restaging the BAM. DNA QC always consumes the NoDup BAM/index, including when coverage is skipped.
- Canonical-only `@SQ` dictionaries are emitted when every retained mate reference permits it. If a canonical record refers to a mate on a noncanonical contig, unused noncanonical `@SQ` lines are retained to keep RNEXT valid; alignment records and coverage signal are still canonical-only.
- nf-core FastQC and the combined Samtools QC process are sidecar QC readers only; they do not alter downstream TrESFlow outputs.
- nf-core `gatk4/markduplicates` replaces the previous local GATK invocation but preserves `--BARCODE_TAG CB`, `--REMOVE_DUPLICATES false`, index creation, and the historical TrESFlow output names via a normalization adapter.
- nf-core `deeptools/bamcoverage` writes the historical `<split>_NoDup.bw` filename directly. Zero-mapped detection is performed inside `SPLIT_DUPLICATES_DNA`; the original NoDup BAM/BAI reaches bamCoverage only when mapped reads are present, while the historical warning TSV is published otherwise.
- `TRES_REPORT_HTML` renders `TrES_Stats/tres_report.html` and `TrES_Stats/tres_report_metrics.json` from existing TrES stats, STARsolo summaries, samtools outputs, and GATK duplicate metrics. FastQC and samtools QC publish under `TrES_Stats/qc/`, with MultiQC nested at `TrES_Stats/qc/multiqc/`.
- The active core runtime lives under [`scripts/core_runtime/`](/mnt/dataFast/ahrmad/tresflowdir/TrESFlow/scripts/core_runtime).
