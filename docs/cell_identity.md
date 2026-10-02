# Cell identity: full-cell-v1

From splitting onward, `CB` and `XI` contain the same complete cell identifier:

```text
<sample_id>_<group_name>_<sb_oligo_index>_<L1L2L3>
```

`sample_id` is the key under `samples`, and `group_name` is the group key.
L1/L2/L3 retain their existing corrected, concatenated order. The positive
oligo index has at least two decimal digits (`1` → `01`, `12` → `12`,
`100` → `100`). Library name, modality, DNA mark, SB nucleotide sequence and
UMI are excluded. `SB` retains the corrected chemistry-specific sequence;
RNA `UM`/`UR` retain the extracted UMI, and STAR `UB` retains the corrected UMI.

The versioned [lookup TSV](../assets/sb_oligo_lookup.v1.tsv) is the runtime
source of truth. It records biological partitions shared across chemistries,
using the mapping supplied from the corrected `BarcodesRefSimple(2).xlsx`.
Excel and Excel-reading packages are not runtime dependencies. RNA always
uses `rna_and_single_dna_sb`, as does single-tag DNA. Dual-tag DNA uses
`dual_dna_sb`. Sequences are exact table entries; dual sequences are never
computed by truncation or reverse complementation.

For example, index `07` resolves to RNA/single-DNA `GCTA` and dual-DNA `GAC`.
Index `10` resolves to `TGCA` and `GCA`. With equal sample, group and ligations,
these different chemistry sequences produce equal full identifiers. All DNA
marks for a cell share its identifier.

## Samplesheet input

The preferred field is shared by all modalities defined by that sample:

```yaml
samples:
  VTD11_VTR12:
    groups:
      B:
        sb_oligo_indices: ["01", "02", 7]
        mark_barcodes:
          H3K27ac: GCCTCTAT
      K422_B:
        sb_oligo_indices: ["12"]
        mark_barcodes:
          H3K27ac: GCCTCTAT
    # Existing rna/dna read blocks follow here.
```

DNA participation still requires group-level `mark_barcodes`. To define
RNA-only or DNA-only groups within a multimodal sample, use the respective
`rna_sb_barcodes` or `dna_sb_barcodes` sequence field. Existing sequence input
is supported and resolved through the lookup. RNA selects `rna_sb_barcodes`
when present, otherwise `sb_barcodes`; single-tag DNA selects
`dna_sb_barcodes`, otherwise `sb_barcodes`. Legacy dual-tag input requires
explicit three-base `dna_sb_barcodes`, together with the RNA four-base input
when RNA participates. Only the new index field expands chemistries
implicitly. The original Tonsil1 samplesheet therefore remains valid.

The parser rejects:

- Indices combined with any sequence-input field in one group.
- Booleans, floats, nondecimal strings, whitespace, zero/negative indices,
  empty lists, unknown indices and duplicates after normalization, including
  `[1, "01"]`.
- Unknown sequences and sequences of the wrong chemistry length.
- Ambiguous lookup rows: missing/extra/empty columns, duplicate normalized
  indices or sequences in either column, lowercase or non-ACGT sequences,
  and four-base/three-base length violations. The entire table is checked,
  including unused rows.
- Index or sequence collisions across groups in one sample, including
  collisions between separate modality-only groups.
- Different RNA/DNA resolved index sets in a shared group. Order may differ.
- Names that would collapse distinct sample/group namespaces or output
  stems. Names are retained exactly and may contain underscores, dots and
  hyphens, including `K422_A`/`K422_B`; whitespace, path/shell characters and
  reserved `Unknown` names are rejected.

Additional positive indices require only valid unique TSV rows. There is no
12- or 99-index limit. `--sb_oligo_lookup /path/to/lookup.tsv` selects another
complete table; its chemistry constraints still apply.

## Derived files and task hashing

`pipeline_info/derived_contract/` contains an exact copy of
`sb_oligo_lookup.v1.tsv` and `cell_identity_version.txt`, including the source
lookup SHA-256. Each modality's SB map has these tab-separated columns:

```text
sample  sb_group  sb_bc  oligo_index  modality  chemistry  input_source  sb_injected_base
```

The first three columns retain their original names and meaning for QC and
report readers. `chemistry` is `rna`, `dna_single` or `dna_dual`;
`input_source` is the original samplesheet field path. `sb_injected_base` is
`-` for ordinary upstream tags or the exact configured A/C/G/T base from
`first_pass_withBC_<base>`. Only this explicit format permits removal of a
leading base, and only from that chemistry's expected raw SB tag. Splitting
requires CB to start with the complete upstream SB, as produced by Tag_Lig3;
there is no blind drop-first lookup or guessed prefix.

`dna_mo_map.tsv` retains `sample`, `sb_group`, `mark`, `mo_bc`.
DNA modality whitelists and `input_fastq_provenance.tsv` retain their contracts.
Metadata carries `cell_id_version`, `sb_lookup_sha256`, complete resolved
`cell_identity_records`, and `split_targets`, an explicit output stem →
sample/group/mark mapping. Group and mark names are never recovered by
splitting on underscores.

The derived maps are Nextflow `path` inputs to tagging, splitting and QC.
The entire lookup's content digest and resolved identities are `val` metadata
through both active graphs, including STAR and nf-core GATK. A lookup or
identity change therefore changes hashes even at an unchanged pathname and
with unchanged size/mtime. Identical derived files are preserved byte-for-byte
without rewriting their timestamps. Cleanup passes the work directory as a
stable string rather than a filesystem object whose timestamp changes during
a run. The resume regression checks cached
unchanged tasks and rerun changed tasks using actual Nextflow traces.

## RNA and DNA processing

RNA FqToSAM carries full CB/XI and separate UM/UR attributes. It never
concatenates names and UMIs in CR or estimates a UMI offset from ID length.
STAR 2.7.11b is configured with `--soloCBtype String`,
`--soloInputSAMattrBarcodeSeq XI UM`, and `--soloCBwhitelist None`.
`--soloCBlen 24` remains a valid parameter-validation placeholder; String
mode reads the two attributes independently. Counting, UMI deduplication and
cell calling distinguish oligo indices before any aggregation. Production
alignment, gene, UMI and EmptyDrops settings remain intact.

Raw/filtered barcode TSVs, matrices and per-cell statistics consequently use
full IDs directly, retaining STAR's matrix column order. STAR preserves input
SAM attributes and emits CB alongside UB, so a streaming validation step
removes identical repeated CB tags and requires CB = XI on every alignment
record. STAR's additional `CB:Z:-` for rejected UMIs is a placeholder;
normalization retains the canonical input CB/XI and leaves UB unchanged,
including `UB:Z:-`. Genuine conflicting non-placeholder identities still
fail. N-containing and homopolymer UMIs retain STAR's existing exclusion from
counts, and normalization preserves aligned reads outside or ambiguous between
genes. The final RNA called-cell filter and retention audit both select XI.

FqToSAM excludes exact `NoMatch` values on CB/SB/MO barcode attributes on either
mate. A permitted sample/group name containing `NoMatch`, including
`NoMatchControl`, is retained; names are not searched for failure substrings.

Picard 3.4.0's `BARCODE_TAG` validates nucleotide UMIs and rejects namespaced
strings. The active nf-core GATK module therefore uses
`--READ_ONE_BARCODE_TAG CB --READ_TWO_BARCODE_TAG SB`. Both mates already have
full CB and corrected SB before alignment/duplicate marking. Picard hashes
string barcode attributes; the unique short chemistry SB provides a second
key, keeping different indices separate even when their full CB hashes
collide. A real GATK fixture exercises such a collision. Real GATK 4.6.2.0 regression data
prove that equal coordinates/ligations with different indices stay separate,
while within-cell genomic duplicates can span physical units and optical
duplicates remain unit-local. DNA `RG`/`PU` retain the AVITI physical unit,
all units share the logical `LB`, and aligned/MarkedDup/NoDup BAMs retain full
CB/XI. Read names and existing filtering/coverage behavior remain intact.

## Regeneration and bounded validation

Older split FASTQs, STAR matrices/cell calls and duplicate-marked/NoDup BAMs
were computed with ligation-only identity. Relabeling these aggregated outputs
cannot undo merged counts or duplicate decisions. Regenerate from before the
split stage, using an upstream pair that still has the complete corrected SB
and Tag_Lig3's `SB + L1L2L3` CB, or from raw FASTQs. Existing references and
raw inputs are reusable. Old tag/UMI/ligation/trim intermediates are
biologically reusable if complete and present, but the normal pipeline starts
from raw input and this change conservatively invalidates their task hashes;
it does not expose an import interface for these intermediates.

With default `--cleanup_work true`, producer-owned intermediate FASTQs are
removed only after their consumers finish, and successful Nextflow work is
cleaned. Published old split FASTQs have already shortened CB and cannot be
passed directly to the new splitter. Published tag-record QC sidecars cannot
replace FASTQs. Keep `--cleanup_work false` during development when testing
resume or inspecting intermediate tags; use a fresh output/work directory
for regeneration and never resume an old production output as a relabeling
exercise.

Run the generated-fixture evidence in the configured environment:

```bash
python -m pytest -q tests/test_full_cell_identity.py tests/test_aviti_optical_duplicates.py tests/test_cell_identity_resume.py
nf-test test tests/default.nf.test tests/multi_fastq.nf.test tests/canonical_bam_filter.nf.test --profile test --ci
```

The STAR fixture changes EmptyDrops to `TopCells 2` only in its temporary
script copy, checks two cells with the same gene/UMI/ligations, and runs the
production BAM filter/audit. No private data, large references or machine
paths are added to its fixtures.

For a separate real-data validation, create a bounded synchronized raw subset
with [make_cell_identity_subset.py](../bin/make_cell_identity_subset.py), then
launch the ordinary pipeline with explicit separate output/work directories.
The subset utility requires PyYAML only in the development environment,
reads at most the requested pairs per technical read set, preserves role
synchronization and refuses an existing destination. Small real subsets may
not meet production cell-calling thresholds; the real synthetic counting
fixture provides the required cell-separation evidence without changing those
thresholds.
