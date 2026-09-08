"""Focused regression tests for the hierarchical samplesheet contract."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests" / "samplesheets" / "group_specific_modalities.yaml"


def nextflow_command():
    java = shutil.which("java")
    jars = []
    if os.environ.get("CONDA_PREFIX"):
        jars.extend(
            Path(os.environ["CONDA_PREFIX"]).glob(
                "share/nextflow/dist/*/nextflow-*-one.jar"
            )
        )
    jars.extend(Path.home().glob(".nextflow/framework/*/nextflow-*-one.jar"))
    if java and jars:
        return [java, "-jar", str(sorted(jars)[-1])]
    nextflow = shutil.which("nextflow")
    if nextflow:
        return [nextflow]
    pytest.skip("Nextflow is not installed")


def test_group_specific_modalities_derive_independent_maps():
    """Shared sample FASTQs yield modality-specific group maps and MO union."""
    nextflow = nextflow_command()

    with tempfile.TemporaryDirectory(prefix="tresflow_contract_") as tmp:
        outdir = Path(tmp) / "output"
        nxf_home = Path(tmp) / "nxf_home"
        env = os.environ.copy()
        env["NXF_HOME"] = str(nxf_home)
        env.setdefault("NXF_OFFLINE", "true")
        result = subprocess.run(
            [
                *nextflow,
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

def test_dna_tagmentation_is_required():
    nextflow = nextflow_command()

    source = FIXTURE.read_text()
    explicit = "      tagmentation: dual\n"
    assert explicit in source
    source = source.replace(explicit, "", 1)

    with tempfile.NamedTemporaryFile(
        "w",
        suffix=".yaml",
        prefix="missing_tagmentation_",
        dir=FIXTURE.parent,
        delete=False,
    ) as handle:
        handle.write(source)
        temporary_sheet = Path(handle.name)

    try:
        with tempfile.TemporaryDirectory(prefix="tresflow_missing_tagmentation_") as tmp:
            outdir = Path(tmp) / "output"
            env = os.environ.copy()
            env["NXF_HOME"] = str(Path(tmp) / "nxf_home")
            env.setdefault("NXF_OFFLINE", "true")

            result = subprocess.run(
                [
                    *nextflow,
                    "run",
                    str(REPO),
                    "-preview",
                    "-ansi-log",
                    "false",
                    "--samplesheet",
                    str(temporary_sheet),
                    "--outdir",
                    str(outdir),
                ],
                cwd=REPO,
                env=env,
                text=True,
                capture_output=True,
                timeout=120,
            )

            combined = result.stdout + result.stderr
            if result.returncode != 0 and "Could not resolve host" in combined:
                pytest.skip("configured nextflow launcher is unavailable offline")

            assert result.returncode != 0
            assert (
                "Missing required field: "
                "samples.shared_reads.dna.tagmentation"
            ) in combined
    finally:
        temporary_sheet.unlink(missing_ok=True)


def test_pipeline_construction_does_not_invoke_host_python(tmp_path):
    nextflow = nextflow_command()

    poison_dir = tmp_path / "poison-bin"
    poison_dir.mkdir()
    marker = tmp_path / "host-python-was-invoked"
    python3 = poison_dir / "python3"
    python3.write_text(
        f"#!/bin/sh\nprintf invoked > {marker}\nexit 97\n", encoding="utf-8"
    )
    python3.chmod(0o755)

    env = os.environ.copy()
    env["PATH"] = f"{poison_dir}{os.pathsep}{env['PATH']}"
    env["PYTHON3_BIN"] = str(python3)
    env["NXF_HOME"] = str(tmp_path / "nxf-home")
    env["NXF_OFFLINE"] = "true"
    result = subprocess.run(
        [
            *nextflow,
            "run",
            str(REPO),
            "-preview",
            "-ansi-log",
            "false",
            "--samplesheet",
            str(REPO / "assets/samplesheet.example.yaml"),
            "--outdir",
            str(tmp_path / "output"),
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )

    combined = result.stdout + result.stderr
    if result.returncode != 0 and "Could not resolve host" in combined:
        pytest.skip("configured nextflow launcher is unavailable offline")
    assert result.returncode == 0, combined
    assert not marker.exists(), "launch-time host python3 was executed"


def test_legacy_runtime_fields_warn_and_do_not_enter_runtime_contract(tmp_path):
    nextflow = nextflow_command()
    env = os.environ.copy()
    env["NXF_HOME"] = str(tmp_path / "nxf-home")
    env["NXF_OFFLINE"] = "true"
    outdir = tmp_path / "output"
    result = subprocess.run(
        [
            *nextflow,
            "run",
            str(REPO),
            "-preview",
            "-ansi-log",
            "false",
            "--samplesheet",
            str(REPO / "tests/samplesheets/explicit_runtime_tmpdir.yaml"),
            "--outdir",
            str(outdir),
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )

    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert (
        "Deprecated samplesheet runtime fields are ignored: "
        "runtime.env_prefix, runtime.tmpdir"
    ) in combined
    assert (outdir / "pipeline_info/runtime_contract.tsv").read_text() == (
        "tool\tconfigured_path\texists\tcurrently_used\n"
    )
