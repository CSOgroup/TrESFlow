#!/usr/bin/env python3
"""Compare isolated Phase 4 outputs with the unchanged RNA golden contract."""

from __future__ import annotations

import argparse
import difflib
import importlib.util
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
NORMALIZER_PATH = REPO / "tests/regression/normalize_outputs.py"


def load_normalizer():
    spec = importlib.util.spec_from_file_location("phase4_normalizer", NORMALIZER_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def is_family_path(path: str) -> bool:
    return (
        (path.startswith("rna_align/") and path != "rna_align/versions.yml")
        or path.endswith(".rna_filter_retention.tsv")
    )


def family_contract(contract: dict[str, Any]) -> dict[str, Any]:
    selected: dict[str, Any] = {}
    selected["files"] = sorted(
        path for path in contract.get("files", []) if is_family_path(path)
    )
    selected["directories"] = sorted(
        path
        for path in contract.get("directories", [])
        if path == "rna_align" or path.startswith("rna_align/")
    )
    for section in (
        "tables",
        "text",
        "star_matrices",
        "star_vectors",
        "star_summaries",
        "bams",
        "bigwigs",
    ):
        selected[section] = {
            path: value
            for path, value in contract.get(section, {}).items()
            if is_family_path(path)
        }
    return selected


def main(args: argparse.Namespace) -> None:
    normalizer = load_normalizer()
    expected_payload = json.loads(args.expected.read_text(encoding="utf-8"))
    expected = family_contract(
        normalizer.canonicalize_contract_runtime_metadata(expected_payload["contract"])
    )
    actual = family_contract(
        normalizer.capture_contract(args.root, "rna_only", args.samtools)
    )

    if len(actual["bams"]) != 1:
        raise SystemExit("Expected exactly one filtered RNA BAM")
    if len(actual["bigwigs"]) != 3:
        raise SystemExit("Expected exactly three RNA bigWigs")
    if len(actual["star_matrices"]) != 3:
        raise SystemExit("Expected exactly three STARsolo matrices")

    if actual == expected:
        print(f"PASS: isolated RNA family matches {args.expected}")
        return

    expected_lines = json.dumps(expected, indent=2, sort_keys=True).splitlines()
    actual_lines = json.dumps(actual, indent=2, sort_keys=True).splitlines()
    print(
        "\n".join(
            difflib.unified_diff(
                expected_lines,
                actual_lines,
                fromfile=str(args.expected),
                tofile=str(args.root),
                n=3,
            )
        )[:100_000]
    )
    raise SystemExit("Isolated Phase 4 semantic regression mismatch")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--expected", required=True, type=Path)
    parser.add_argument("--samtools", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
