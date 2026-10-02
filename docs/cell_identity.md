# Cell identity: full-cell-v2

From splitting onward, `CB` and `XI` contain the same complete cell identifier:

```text
<sample_id>_<group_name>_<sb_index>_<L1L2L3>
```

`sample_id` is the key under `samples`, and `group_name` is the group key.
L1/L2/L3 retain their existing corrected, concatenated order. The positive
`sb_index` has at least two decimal digits (`1` → `01`, `12` → `12`,
`100` → `100`). Library name, modality, DNA mark, SB nucleotide sequence and
UMI are excluded. `SB` retains the corrected chemistry-specific sequence;
RNA `UM`/`UR` retain the extracted UMI, and STAR `UB` retains the corrected UMI.

The versioned [lookup TSV](../assets/sb_oligo_lookup.v1.tsv) is the runtime
source of truth for physical oligos. It records corresponding physical oligos across chemistries,
using the mapping supplied from the corrected `BarcodesRefSimple(2).xlsx`.
Excel and Excel-reading packages are not runtime dependencies. RNA always
uses `rna_and_single_dna_sb`, as does single-tag DNA. Dual-tag DNA uses
`dual_dna_sb`. Sequences are exact table entries; dual sequences are never
computed by truncation or reverse complementation.

For example, index `07` resolves to RNA/single-DNA `GCTA` and dual-DNA `GAC`.
Index `10` resolves to `TGCA` and `GCA`. With equal sample, group and ligations,
these different chemistry sequences produce equal full identifiers. All DNA
marks for a cell share its identifier. Without pairing, `sb_index` is the
actual physical `oligo_index`, so existing valid samplesheets keep their IDs.
Indices 13–16 resolve to `GATC`, `TGAC`, `CATG`, `TACG` for RNA/single DNA.
Their dual-DNA column is explicitly unavailable (`-`). Unavailable entries
are excluded from reverse lookup and sequence uniqueness checks; all real
sequences and normalized indices are validated, including unused rows.
Selecting an unavailable physical dual oligo fails validation. It is never
inferred, truncated or substituted. Mappings 03=`GACT` and 11=`CTGA` remain intact.

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
`rna_sb_oligo_indices` / `dna_sb_oligo_indices`, or their existing
`rna_sb_barcodes` / `dna_sb_barcodes` sequence fields. Existing sequence input
is supported and resolved through the lookup. RNA selects `rna_sb_barcodes`
when present, otherwise `sb_barcodes`; single-tag DNA selects
`dna_sb_barcodes`, otherwise `sb_barcodes`. Legacy dual-tag input requires
explicit three-base `dna_sb_barcodes`, together with the RNA four-base input
when RNA participates. Only index selectors expand chemistries
implicitly. The original Tonsil1 samplesheet therefore remains valid.

Independent modality selections can differ. The parser warns and continues,
retaining each actual index in its identifier. Matching sets in different
orders are equivalent; list order never establishes correspondence.

```yaml
B:
  rna_sb_oligo_indices: ["07", "13"]
  dna_sb_oligo_indices: ["10", "07"]
  mark_barcodes:
    H3K27ac: GCCTCTAT
```

Shared `sb_oligo_indices` is exclusive with all other barcode selectors.
Modality index lists may be mixed with the *other* modality's index or sequence
input, but conflict with their own modality's sequence field or shared
`sb_barcodes` when it applies to that modality. Dual DNA indices may coexist
with shared `sb_barcodes` used for RNA. Existing sequence-only precedence remains intact.

### Exceptional explicit pairing

Use `sb_oligo_pairings` only when the samplesheet must explicitly correct an
experimenter mismatch between physical oligos representing one biological
partition. Standard corresponding selections use the normal interface.

```yaml
B:
  sb_oligo_pairings:
    - sb_index: "07"
      rna_oligo_index: "07"
      dna_oligo_index: "10"
    - sb_index: "13"
      rna_sb_barcode: GATC
      dna_sb_barcode: GAC
  mark_barcodes:
    H3K27ac: GCCTCTAT
```

For dual DNA, the first entry selects RNA `GCTA` (physical 07) and DNA `GCA`
(physical 10), both labelled `07`. The second selects RNA physical 13 and
DNA physical 07, both labelled `13`. This succeeds even though physical 13
has no dual sequence: the selected dual oligo is 07. `sb_index` is an identity
label and is never used to select chemistry sequences. Any strictly positive
decimal label is permitted, including one absent from the physical lookup.

Each entry requires `sb_index` and at least one modality. Each present modality
requires exactly one of `rna_oligo_index` / `rna_sb_barcode`, or
`dna_oligo_index` / `dna_sb_barcode`. Sequence selectors are singular exact
chemistry-specific lookup sequences. Omission means a modality-specific
partition; it does not create an inferred counterpart. Complementary entries
may share a label, but each shared identity permits at most one physical
oligo per modality. All selectors require the corresponding sample reads
block. Pairing is exclusive with every other group barcode-selection field.

Pairing applies only within the same sample and group. It changes CB/XI and
matrix/per-cell identity labels while preserving physical SB, molecular UMIs,
DNA marks and AVITI read groups. Physical `oligo_index` remains independent
of logical `sb_index` in derived records and audits.

After the complete samplesheet validates, every nontrivial remapping emits a
clearly delimited Nextflow warning before runtime/reference preflight or any
processing. It names the sample/group, actual indices and sequences, selected
DNA chemistry and shared `sb_index`, and states that the samplesheet declares
one biological partition. Independent mismatches warn that unmatched full
identifiers differ. Cross-modality reuse across different groups warns about
distinct namespaces when routing remains unambiguous. Warnings are emitted
again on resumed runs, require no confirmation and do not pause processing.

The parser rejects:

- Conflicting selectors, missing/extra pairing fields, and reused physical
  barcodes or multiple physical oligos per shared identity within a modality.
- Booleans, floats, nondecimal strings, whitespace, zero/negative indices,
  empty lists, unknown physical indices and duplicates after normalization.
- Unknown sequences, wrong chemistry lengths, and unavailable selected dual
  oligos. The logical label alone never requires a physical lookup row.
- Malformed lookup rows, duplicate normalized indices, duplicate real
  sequences, lowercase/non-ACGT sequences and length violations. Only the
  explicit `-` in the dual column is an unavailable sentinel.
- Physical barcode collisions across groups within one sample and modality.
- Names that collapse distinct sample/group namespaces or output stems.
  Names may retain underscores, dots and hyphens; namespaces are explicit
  metadata and are never inferred by splitting underscores.

Additional positive physical indices require valid unique TSV rows. There is
no 12- or 99-index limit. `--sb_oligo_lookup /path/to/lookup.tsv` selects another
complete table. The [samplesheet schema](../assets/samplesheet.schema.json)
covers selector structure for editor/tool validation; the parser additionally
checks normalized uniqueness, chemistry availability and routing.

## Derived files and task hashing

`pipeline_info/derived_contract/` contains an exact copy of
`sb_oligo_lookup.v1.tsv` and `cell_identity_version.txt`, including the source
lookup SHA-256. Each modality's SB map has these tab-separated columns:

```text
sample  sb_group  sb_bc  oligo_index  modality  chemistry  input_source  sb_injected_base  sb_index
```

`oligo_index` retains the actual physical lookup index; the appended
`sb_index` supplies the logical identity used by both splitters. Older maps
without that appended field default to their physical index.
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
`pipeline_info/sb_identity_warnings.txt` saves readable warnings (or explicitly
records no warnings). `pipeline_info/sb_physical_to_logical.tsv` audits all
resolved modality mappings with the same columns as the SB maps. These files
are regenerated after successful validation on every launch, including resume.

Metadata carries `cell_id_version`, `sb_lookup_sha256`, `sb_mapping_sha256`, complete resolved
`cell_identity_records`, and `split_targets`, an explicit output stem →
sample/group/mark mapping. Group and mark names are never recovered by
splitting on underscores.

The derived maps are Nextflow `path` inputs to tagging, splitting and QC.
The entire lookup's content digest and resolved identities are `val` metadata
through both active graphs, including STAR and nf-core GATK. A lookup or
identity change, including a change only to `sb_index` or pairing, therefore changes hashes even at an unchanged pathname and
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
cell calling distinguish logical SB indices before any aggregation. Production
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
python -m pytest -q tests/test_full_cell_identity.py tests/test_sb_index_pairings.py tests/test_aviti_optical_duplicates.py tests/test_cell_identity_resume.py
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
