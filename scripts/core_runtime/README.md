# Core Runtime Scripts

This directory contains the active repo-owned runtime for the validated core workflow.

Purpose:

- keep the implemented core pipeline runnable without runtime dependence on `upstream/source_scripts/`
- preserve the validated RNA and DNA behavior in repo-owned code
- keep `upstream/source_scripts/` only as provenance

The current core runtime set is:

- `Tag.codon`
- `Tag_UMI.codon`
- `Tag_Lig3.codon`
- `Split_ReadsV2.codon`
- `FqToSAM.codon`
- `RNA_STARSOLO_ALIGN.sh`
- `NormalizeRnaBamTags.py`
- `RNA_FILTERED_BAM.sh`
- `RNA_COVERAGE.sh`
- `AlignDNA.sh`
- `FilterCanonicalBam.sh`
- `utils.codon`

Ownership rules:

- edits to files in this directory are pipeline changes and should be reviewed like any other repo code
- keep behavior aligned with the validated workflow unless a deliberate contract change is documented
- prefer targeted readability improvements over large rewrites of biological logic

Splitting and RNA SAM transport implement [full-cell-v2](../../docs/cell_identity.md).
SB identity maps retain a normalized physical `oligo_index` and append logical
`sb_index` (defaulting to the physical index in older maps). Splitters construct
CB/XI from `sb_index` and retain the actual corrected SB sequence. STAR uses String
barcodes from separate XI/UM attributes; it never slices a namespaced ID to
find a UMI. BAM tag normalization checks every alignment class.
