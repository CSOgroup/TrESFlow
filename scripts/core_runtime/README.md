# Core Runtime Scripts

This directory contains the **repo-owned runtime copies** of the currently
validated core workflow scripts.

Purpose:

- remove runtime dependence on `upstream/source_scripts/` for the implemented core pipeline
- preserve the validated behavior of the current RNA and DNA branches
- keep the upstream source tree available only for provenance and comparison

The current core runtime set is:

- `Tag.codon`
- `Tag_UMI.codon`
- `Tag_Lig3.codon`
- `Split_ReadsV2.codon`
- `FqToSAM.codon`
- `AlignRNA.sh`
- `AlignDNA.sh`
- `utils.codon`

Ownership rules:

- edits to files in this directory are pipeline changes and should be reviewed like any other repo code
- keep behavior aligned with the validated workflow unless a deliberate contract change is documented
- do not point the core Nextflow modules back at `upstream/source_scripts/` except for temporary debugging

The optional downstream `sc_process.py` path is intentionally separate and may
still use the upstream source tree until that optional component is refactored
independently.
