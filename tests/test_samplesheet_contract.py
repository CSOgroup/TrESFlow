"""Focused regression tests for the hierarchical samplesheet contract."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests" / "samplesheets" / "group_specific_modalities.yaml"


def test_group_specific_modalities_derive_independent_maps():
    """Shared sample FASTQs yield modality-specific group maps and MO union."""
    nextflow = shutil.which("nextflow")
    if not nextflow:
        pytest.skip("nextflow is not installed")

    with tempfile.TemporaryDirectory(prefix="tresflow_contract_") as tmp:
        outdir = Path(tmp) / "output"
        nxf_home = Path(tmp) / "nxf_home"
        env = os.environ.copy()
        env["NXF_HOME"] = str(nxf_home)
        env.setdefault("NXF_OFFLINE", "true")
        result = subprocess.run(
            [
                nextflow,
                "run",
                str(REPO),
                "-preview",
                "-ansi-log",
                "false",
                "--samplesheet",
                str(FIXTURE),
                "--outdir",
                str(outdir),
            ],
            cwd=REPO,
            env=env,
            text=True,
            capture_output=True,
            timeout=120,
        )
        if result.returncode != 0 and "Could not resolve host" in (result.stdout + result.stderr):
            pytest.skip("configured nextflow launcher is unavailable offline")
        assert result.returncode == 0, result.stdout + result.stderr

        contract = outdir / "pipeline_info" / "derived_contract"
        rna_rows = (contract / "rna_sb_group_map.tsv").read_text().splitlines()[1:]
        dna_rows = (contract / "dna_sb_group_map.tsv").read_text().splitlines()[1:]
        mo_rows = (contract / "dna_mo_map.tsv").read_text().splitlines()[1:]
        whitelist = (contract / "dna_modality_whitelists" / "shared_reads.txt").read_text().split()

        assert {row.split("\t")[1] for row in rna_rows} == {"RNA_only"}
        assert {row.split("\t")[1] for row in dna_rows} == {"DNA_only", "DNA_only_2"}
        assert {tuple(row.split("\t")[1:]) for row in mo_rows} == {
            ("DNA_only", "MarkA", "AGGCTATA"),
            ("DNA_only", "MarkB", "GCCTCTAT"),
            ("DNA_only_2", "MarkB", "AGGCTATA"),
        }
        assert set(whitelist) == {"AGGCTATA", "GCCTCTAT"}
