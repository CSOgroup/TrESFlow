#!/usr/bin/env python3
"""Compare isolated Phase 6 QC outputs with the unchanged v1.1.1 contracts."""

from __future__ import annotations

import argparse
import csv
import difflib
import importlib.util
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
NORMALIZER_PATH = REPO / "tests/regression/normalize_outputs.py"
EXPECTED_CONTAINERS = {
    "FASTQC": (
        "quay.io/biocontainers/fastqc@"
        "sha256:e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900"
    ),
    "SAMTOOLS_BAM_QC": (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    "TRES_REPORT_HTML": (
        "docker.io/library/python:3.12.13-bookworm@"
        "sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12"
    ),
    "MULTIQC": (
        "quay.io/biocontainers/multiqc@"
        "sha256:dfd9fde2c48b896b884e79a71ddc16c72c97a0ee5c5c8e45aaba50f55d07d263"
    ),
}
REPORT_FILES = {
    "TrES_Stats/tres_report.html",
    "TrES_Stats/read_retention.tsv",
    "TrES_Stats/qc_metrics.tsv",
    "TrES_Stats/barcode_composition.tsv",
    "TrES_Stats/library_complexity.tsv",
}
FORBIDDEN_TASK_COUPLING = (
    "/home/annan/micromamba/envs/tres",
    "${projectDir}",
    "runtime.env_prefix",
    "PYTHON3_BIN",
    "SAMTOOLS_BIN",
)


def load_normalizer():
    spec = importlib.util.spec_from_file_location("phase6_normalizer", NORMALIZER_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def is_qc_path(path: str) -> bool:
    return path.startswith("TrES_Stats/qc/") or path in REPORT_FILES


def is_optional_plot_file(path: str) -> bool:
    return path.startswith("TrES_Stats/qc/multiqc/multiqc_report_plots/")


def qc_contract(contract: dict[str, Any]) -> dict[str, Any]:
    selected: dict[str, Any] = {
        "directories": sorted(
            path for path in contract.get("directories", []) if is_qc_path(path)
        ),
        "files": sorted(
            path
            for path in contract.get("files", [])
            if is_qc_path(path) and not is_optional_plot_file(path)
        )
    }
    for section in ("tables", "text", "json", "gzip_text", "html", "fastqc", "multiqc"):
        value = contract.get(section, {})
        if section == "multiqc":
            selected[section] = value
        else:
            selected[section] = {
                path: item
                for path, item in value.items()
                if is_qc_path(path) and not is_optional_plot_file(path)
            }
    return selected


def validate_optional_plots(args: argparse.Namespace) -> None:
    plots_root = args.root / "TrES_Stats/qc/multiqc/multiqc_report_plots"
    if not plots_root.is_dir():
        raise AssertionError("MultiQC plot publication directory is missing")
    observed: dict[str, set[str]] = {}
    for extension in ("pdf", "png", "svg"):
        directory = plots_root / extension
        if not directory.is_dir():
            raise AssertionError(f"MultiQC {extension} plot directory is missing")
        files = list(directory.glob(f"*.{extension}"))
        if any(path.stat().st_size == 0 for path in files):
            raise AssertionError(f"MultiQC emitted an empty {extension} plot")
        observed[extension] = {path.stem for path in files}
    nonempty = [names for names in observed.values() if names]
    if nonempty and any(names != nonempty[0] for names in nonempty[1:]):
        raise AssertionError(f"MultiQC plot formats do not cover the same plots: {observed}")


def process_name(row: dict[str, str]) -> str:
    return row["name"].split(" (", 1)[0].rsplit(":", 1)[-1]


def validate_trace(args: argparse.Namespace) -> None:
    expected_counts = {
        "FASTQC": 1,
        "SAMTOOLS_BAM_QC": 1 if args.scenario == "rna_only" else 3,
        "TRES_REPORT_HTML": 1,
        "MULTIQC": 1,
    }
    with args.trace.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    observed = {name: 0 for name in expected_counts}
    for row in rows:
        name = process_name(row)
        if name not in observed:
            raise AssertionError(f"Unexpected Phase 6 task in trace: {row['name']}")
        observed[name] += 1
        if row["status"] != "COMPLETED" or row["exit"] != "0":
            raise AssertionError(f"Task did not complete: {row}")
        if row["container"] != EXPECTED_CONTAINERS[name]:
            raise AssertionError(
                f"Unexpected container declaration for {row['name']}: {row['container']}"
            )

        workdir = Path(row["workdir"])
        script = (workdir / ".command.sh").read_text(encoding="utf-8")
        wrapper = (workdir / ".command.run").read_text(encoding="utf-8")
        for token in FORBIDDEN_TASK_COUPLING:
            if token in script:
                raise AssertionError(f"Forbidden host coupling {token!r}: {row['name']}")
        if name in {"FASTQC", "SAMTOOLS_BAM_QC", "MULTIQC"}:
            if 'export TMPDIR="$PWD/.tmp"' not in script:
                raise AssertionError(f"Task-local TMPDIR is missing: {row['name']}")
        if args.engine == "docker":
            if "docker run " not in wrapper or EXPECTED_CONTAINERS[name] not in wrapper:
                raise AssertionError(f"Task was not executed by Docker: {row['name']}")
        else:
            if "# conda environment" not in wrapper or "micromamba activate " not in wrapper:
                raise AssertionError(f"Task did not activate Conda: {row['name']}")
            for executor in ("docker run ", "apptainer exec ", "singularity exec "):
                if executor in wrapper:
                    raise AssertionError(
                        f"Conda task unexpectedly used {executor.strip()}: {row['name']}"
                    )
    if observed != expected_counts:
        raise AssertionError(f"Unexpected Phase 6 task counts: {observed}")


def compare_contract(args: argparse.Namespace) -> None:
    normalizer = load_normalizer()
    payload = json.loads((args.expected_dir / f"{args.scenario}.json").read_text())
    expected = qc_contract(
        normalizer.canonicalize_contract_runtime_metadata(payload["contract"])
    )
    actual = qc_contract(
        normalizer.capture_contract(args.root, args.scenario, args.samtools)
    )
    if actual != expected:
        expected_lines = json.dumps(expected, indent=2, sort_keys=True).splitlines()
        actual_lines = json.dumps(actual, indent=2, sort_keys=True).splitlines()
        print(
            "\n".join(
                difflib.unified_diff(
                    expected_lines,
                    actual_lines,
                    fromfile=str(args.expected_dir / f"{args.scenario}.json"),
                    tofile=str(args.root),
                    n=3,
                )
            )[:100_000]
        )
        raise AssertionError(
            f"Isolated Phase 6 {args.scenario} semantic regression mismatch"
        )
    print(
        f"PASS: Phase 6 {args.engine} {args.scenario} QC matches the v1.1.1 contract"
    )


def main(args: argparse.Namespace) -> None:
    validate_trace(args)
    validate_optional_plots(args)
    compare_contract(args)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--engine", required=True, choices=("docker", "conda"))
    parser.add_argument(
        "--scenario", required=True, choices=("rna_only", "dna_single", "dna_dual")
    )
    parser.add_argument(
        "--expected-dir",
        type=Path,
        default=REPO / "tests/regression/contracts/v1.1.1",
    )
    parser.add_argument("--samtools", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
