#!/usr/bin/env python3
"""Run real Phase 0 scenarios against v1.1.1 and the current worktree."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path


SCENARIOS = ("rna_only", "dna_single", "dna_dual")
BASELINE_COMMIT = "40a383e80d945618952f8e5bfddb73ed3dc63af6"


def run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    print("+ " + " ".join(command), flush=True)
    result = subprocess.run(command, cwd=cwd, env=env)
    if result.returncode:
        raise SystemExit(result.returncode)


def output(command: list[str], *, cwd: Path) -> str:
    return subprocess.run(
        command, cwd=cwd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_asset_manifest(fixture_root: Path) -> None:
    manifest_path = fixture_root / "asset-manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    errors = []
    for entry in manifest["files"]:
        path = fixture_root / entry["path"]
        if not path.is_file():
            errors.append(f"missing: {entry['path']}")
            continue
        if path.stat().st_size != entry["bytes"]:
            errors.append(f"size: {entry['path']}")
        if sha256(path) != entry["sha256"]:
            errors.append(f"sha256: {entry['path']}")
    if errors:
        raise SystemExit("Fixture manifest validation failed:\n" + "\n".join(errors))


def validate_baseline_provenance(expected_dir: Path) -> dict[str, object]:
    provenance_path = expected_dir / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("release") != "v1.1.1" or provenance.get("commit") != BASELINE_COMMIT:
        raise SystemExit(f"Unexpected baseline provenance in {provenance_path}")
    for scenario in SCENARIOS:
        expected_path = expected_dir / f"{scenario}.json"
        if not expected_path.is_file():
            continue
        payload = json.loads(expected_path.read_text(encoding="utf-8"))
        if payload.get("baseline_provenance") != provenance:
            raise SystemExit(f"Baseline provenance mismatch in {expected_path}")
    return provenance


def safe_extract(archive_path: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive_path) as archive:
        for member in archive.getmembers():
            target = (destination / member.name).resolve()
            if target != root and root not in target.parents:
                raise SystemExit(f"Unsafe path in git archive: {member.name}")
        archive.extractall(destination, filter="data")


def pipeline_command(
    args: argparse.Namespace,
    project: Path,
    fixture_root: Path,
    scenario: str,
    output_dir: Path,
    work_dir: Path,
) -> list[str]:
    command = [
        str(args.nextflow),
        "run",
        str(project),
    ]
    if args.engine != "standard":
        command.append(f"-with-{args.engine}")
    return [
        *command,
        "-profile",
        "standard",
        "-c",
        str(args.repo / "tests/regression/nextflow.config"),
        "--samplesheet",
        str(fixture_root / "samplesheets" / f"{scenario}.yaml"),
        "--outdir",
        str(output_dir),
        "--max_cpus",
        "4",
        "--cleanup_work",
        "false",
        "-work-dir",
        str(work_dir),
    ]


def normalizer_command(args: argparse.Namespace) -> list[str]:
    return [
        str(args.python),
        str(args.repo / "tests/regression/normalize_outputs.py"),
        "--samtools",
        str(args.samtools),
    ]


def main(args: argparse.Namespace) -> None:
    args.repo = args.repo.resolve()
    args.workspace = args.workspace.resolve()
    args.expected_dir = args.expected_dir.resolve()
    if args.workspace.exists() and any(args.workspace.iterdir()):
        raise SystemExit(f"Refusing to reuse non-empty regression workspace: {args.workspace}")
    args.workspace.mkdir(parents=True, exist_ok=True)
    validate_baseline_provenance(args.expected_dir)

    baseline_commit = output(["git", "rev-parse", "v1.1.1^{}"], cwd=args.repo)
    if baseline_commit != BASELINE_COMMIT:
        raise SystemExit(f"v1.1.1 resolved to {baseline_commit}, expected {BASELINE_COMMIT}")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", "v1.1.1", "HEAD"], cwd=args.repo
    )
    if ancestor.returncode:
        raise SystemExit("v1.1.1 is not an ancestor of HEAD")

    fixture_root = args.workspace / "fixture"
    run(
        [
            str(args.python),
            str(args.repo / "tests/regression/generate_fixture.py"),
            "--output",
            str(fixture_root),
            "--env-prefix",
            str(args.env_prefix),
            "--star",
            str(args.star),
            "--bwa-mem2",
            str(args.bwa_mem2),
        ],
        cwd=args.repo,
    )
    validate_asset_manifest(fixture_root)

    baseline_project = args.workspace / "source-v1.1.1"
    archive = args.workspace / "v1.1.1.tar"
    run(
        ["git", "archive", "--format=tar", "--output", str(archive), "v1.1.1"],
        cwd=args.repo,
    )
    safe_extract(archive, baseline_project)
    extracted_version = output(
        ["git", "show", "v1.1.1:nextflow.config"], cwd=args.repo
    )
    if "version         = 'v1.1.1'" not in extracted_version:
        raise SystemExit("The extracted baseline does not declare pipeline version v1.1.1")

    run_env = os.environ.copy()
    run_env.update({"NXF_OFFLINE": "true", "NXF_ANSI_LOG": "false"})
    scenarios = SCENARIOS if args.scenario == "all" else (args.scenario,)
    provenance = args.expected_dir / "provenance.json"

    for scenario in scenarios:
        print(f"\n=== {scenario}: v1.1.1 ===", flush=True)
        baseline_output = args.workspace / "outputs" / "v1.1.1" / scenario
        baseline_work = args.workspace / "work" / "v1.1.1" / scenario
        baseline_env = run_env | {"NXF_HOME": str(args.workspace / "nxf-home" / "v1.1.1")}
        run(
            pipeline_command(
                args, baseline_project, fixture_root, scenario, baseline_output, baseline_work
            ),
            cwd=baseline_project,
            env=baseline_env,
        )

        expected = args.expected_dir / f"{scenario}.json"
        if args.capture_baseline:
            run(
                [
                    *normalizer_command(args),
                    "capture",
                    "--root",
                    str(baseline_output),
                    "--scenario",
                    scenario,
                    "--provenance",
                    str(provenance),
                    "--output",
                    str(expected),
                ],
                cwd=args.repo,
            )
        if not expected.is_file():
            raise SystemExit(
                f"Missing golden contract {expected}; capture requires the explicit --capture-baseline flag"
            )
        run(
            [*normalizer_command(args), "compare", "--root", str(baseline_output), "--expected", str(expected)],
            cwd=args.repo,
        )

        print(f"\n=== {scenario}: current worktree ===", flush=True)
        current_output = args.workspace / "outputs" / "current" / scenario
        current_work = args.workspace / "work" / "current" / scenario
        current_env = run_env | {"NXF_HOME": str(args.workspace / "nxf-home" / "current")}
        run(
            pipeline_command(args, args.repo, fixture_root, scenario, current_output, current_work),
            cwd=args.repo,
            env=current_env,
        )
        run(
            [*normalizer_command(args), "compare", "--root", str(current_output), "--expected", str(expected)],
            cwd=args.repo,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--env-prefix", required=True, type=Path)
    parser.add_argument("--nextflow", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--star", required=True, type=Path)
    parser.add_argument("--bwa-mem2", required=True, type=Path)
    parser.add_argument("--samtools", required=True, type=Path)
    parser.add_argument(
        "--engine",
        choices=("standard", "docker", "conda", "apptainer"),
        default="standard",
        help="Nextflow task execution engine; use a separate empty workspace for each engine",
    )
    parser.add_argument("--scenario", choices=(*SCENARIOS, "all"), default="all")
    parser.add_argument(
        "--expected-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "contracts/v1.1.1",
    )
    parser.add_argument(
        "--capture-baseline",
        action="store_true",
        help="Create missing golden files from v1.1.1; existing files are never overwritten",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
