#!/usr/bin/env python3
"""Compare isolated Phase 5 outputs with the unchanged DNA golden contracts."""

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
    "ALIGN_DNA": (
        "community.wave.seqera.io/library/bwa-mem2_samtools@"
        "sha256:ce8cbf5cc21c690c8c2994d9bbb409b9313c47f38c5452a5fe7dec3402eff9c8"
    ),
    "FILTER_CANONICAL_DNA_ALIGNED_BAM": (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    "GATK4_MARKDUPLICATES": (
        "community.wave.seqera.io/library/gatk4_samtools@"
        "sha256:ca703ab322c8cf829f3987275648ab98fec87b4a340b5cbe154a49cb9a4f41a0"
    ),
    "NORMALIZE_DNA_MARKDUPLICATES": (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    "SPLIT_DUPLICATES_DNA": (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    "DEEPTOOLS_BAMCOVERAGE": (
        "community.wave.seqera.io/library/deeptools_samtools@"
        "sha256:9c482aa632f9a30dcd09c38bed01ac89c54b345fd1ec666deb3ec9c595b258c9"
    ),
}
EXPECTED_TASK_COUNTS = {
    "ALIGN_DNA": 2,
    "FILTER_CANONICAL_DNA_ALIGNED_BAM": 2,
    "GATK4_MARKDUPLICATES": 2,
    "NORMALIZE_DNA_MARKDUPLICATES": 2,
    "SPLIT_DUPLICATES_DNA": 3,
    "DEEPTOOLS_BAMCOVERAGE": 2,
}
FORBIDDEN_TASK_COUPLING = (
    "/home/annan/micromamba/envs/tres",
    "${projectDir}",
    "runtime.env_prefix",
    "BWA_MEM2_BIN",
    "GATK_BIN",
    "BAMCOVERAGE_BIN",
)


def load_normalizer():
    spec = importlib.util.spec_from_file_location("phase5_normalizer", NORMALIZER_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def is_family_path(path: str, split_name: str) -> bool:
    return path in {
        f"TrES_Stats/{split_name}.dna_alignment_retention.tsv",
        f"dna_align/{split_name}.DuplicateMetrics.txt",
        f"dna_align/{split_name}_MarkedDup.bam",
        f"dna_align/{split_name}_MarkedDup.bam.bai",
        f"dna_align/{split_name}_NoDup.bam",
        f"dna_align/{split_name}_NoDup.bam.bai",
        f"dna_align/{split_name}_NoDup.bw",
    }


def family_contract(contract: dict[str, Any], split_name: str) -> dict[str, Any]:
    selected: dict[str, Any] = {
        "files": sorted(
            path
            for path in contract.get("files", [])
            if is_family_path(path, split_name)
        )
    }
    for section in ("tables", "text", "bams", "bigwigs"):
        selected[section] = {
            path: value
            for path, value in contract.get(section, {}).items()
            if is_family_path(path, split_name)
        }
    return selected


def assert_bam_invariants(contract: dict[str, Any], sample_id: str) -> None:
    marked_path = next(path for path in contract["bams"] if path.endswith("_MarkedDup.bam"))
    nodup_path = next(path for path in contract["bams"] if path.endswith("_NoDup.bam"))
    marked = contract["bams"][marked_path]
    nodup = contract["bams"][nodup_path]

    expected_rg = (
        "@RG\tID:AVTEST:RUN1:FLOW1:L1\tLB:PHASE0_SYNTHETIC\tPL:ELEMENT"
        f"\tPM:AVITI_500MIO\tPU:AVTEST:RUN1:FLOW1:L1\tSM:{sample_id}"
    )
    for label, bam in (("MarkedDup", marked), ("NoDup", nodup)):
        read_groups = [line for line in bam["header"] if line.startswith("@RG\t")]
        if read_groups != [expected_rg]:
            raise AssertionError(f"{label} did not preserve the lane-only AVITI @RG: {read_groups}")
        for record in bam["records"]:
            tags = {field.split(":", 1)[0] for field in record[11:]}
            missing = {"CB", "RG", "SB", "MO"} - tags
            if missing:
                raise AssertionError(f"{label} record {record[0]} lacks tags {sorted(missing)}")

    if not any(int(record[1]) & 0x400 for record in marked["records"]):
        raise AssertionError("MarkedDup BAM does not exercise duplicate classification")
    if any(int(record[1]) & 0x400 for record in nodup["records"]):
        raise AssertionError("NoDup BAM retains records carrying flag 0x400")


def picard_metrics(lines: list[str]) -> dict[str, str]:
    for index, line in enumerate(lines):
        if not line.startswith("LIBRARY\t"):
            continue
        header = line.split("\t")
        for values_line in lines[index + 1 :]:
            if values_line and not values_line.startswith("#"):
                values = values_line.split("\t")
                values += [""] * (len(header) - len(values))
                return dict(zip(header, values, strict=True))
    raise AssertionError("Picard duplicate metrics table is missing")


def assert_duplicate_metrics(contract: dict[str, Any]) -> None:
    metrics_path = next(
        path for path in contract["text"] if path.endswith(".DuplicateMetrics.txt")
    )
    metrics = picard_metrics(contract["text"][metrics_path])
    required = {
        "READ_PAIRS_EXAMINED",
        "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES",
        "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE",
    }
    if not required.issubset(metrics):
        raise AssertionError(f"Picard metrics omit required fields: {sorted(required - metrics.keys())}")
    if int(metrics["READ_PAIRS_EXAMINED"]) <= 0:
        raise AssertionError("Picard edge fixture has no examined read pairs")
    if int(metrics["READ_PAIR_DUPLICATES"]) <= 0:
        raise AssertionError("Picard edge fixture has no duplicate read pairs")
    if int(metrics["READ_PAIR_OPTICAL_DUPLICATES"]) <= 0:
        raise AssertionError("Picard edge fixture has no optical duplicate read pairs")


def process_name(row: dict[str, str]) -> str:
    return row["name"].split(" (", 1)[0].rsplit(":", 1)[-1]


def validate_trace(args: argparse.Namespace) -> None:
    with args.trace.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    observed = {name: 0 for name in EXPECTED_TASK_COUNTS}
    for row in rows:
        name = process_name(row)
        if name not in observed:
            raise AssertionError(f"Unexpected Phase 5 task in trace: {row['name']}")
        observed[name] += 1
        if row["status"] != "COMPLETED" or row["exit"] != "0":
            raise AssertionError(f"Task did not complete: {row}")
        if row["container"] != EXPECTED_CONTAINERS[name]:
            raise AssertionError(f"Unexpected container declaration for {row['name']}: {row['container']}")

        workdir = Path(row["workdir"])
        script = (workdir / ".command.sh").read_text(encoding="utf-8")
        wrapper = (workdir / ".command.run").read_text(encoding="utf-8")
        for token in FORBIDDEN_TASK_COUPLING:
            if token in script:
                raise AssertionError(f"Forbidden host coupling {token!r}: {row['name']}")
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
                    raise AssertionError(f"Conda task unexpectedly used {executor.strip()}: {row['name']}")

    if observed != EXPECTED_TASK_COUNTS:
        raise AssertionError(f"Unexpected Phase 5 task counts: {observed}")

    coverage_rows = [row for row in rows if process_name(row) == "DEEPTOOLS_BAMCOVERAGE"]
    if len(coverage_rows) != 2:
        raise AssertionError("Expected exactly the two nonempty NoDup BAMs to reach bamCoverage")
    for row in coverage_rows:
        script = (Path(row["workdir"]) / ".command.sh").read_text(encoding="utf-8")
        bam_lines = [line.strip() for line in script.splitlines() if line.strip().startswith("--bam ")]
        if len(bam_lines) != 1 or "_NoDup.bam" not in bam_lines[0]:
            raise AssertionError(f"bamCoverage input is not a NoDup BAM: {row['name']} {bam_lines}")
        if "_MarkedDup.bam" in bam_lines[0] or "phase5_empty" in bam_lines[0]:
            raise AssertionError(f"Invalid bamCoverage source: {row['name']} {bam_lines[0]}")


def compare_scenario(args: argparse.Namespace, scenario: str) -> None:
    sample_id = f"phase0_{scenario}"
    split_name = f"{sample_id}_Normal_H3K27ac"
    expected_path = args.expected_dir / f"{scenario}.json"
    normalizer = load_normalizer()
    payload = json.loads(expected_path.read_text(encoding="utf-8"))
    expected = family_contract(
        normalizer.canonicalize_contract_runtime_metadata(payload["contract"]),
        split_name,
    )
    actual = family_contract(
        normalizer.capture_contract(args.root, scenario, args.samtools),
        split_name,
    )

    if len(actual["bams"]) != 2 or len(actual["bigwigs"]) != 1:
        raise AssertionError(f"Incomplete {scenario} DNA family capture")
    if not next(iter(actual["bigwigs"])).endswith("_NoDup.bw"):
        raise AssertionError(f"{scenario} coverage output is not named from NoDup")
    assert_bam_invariants(actual, sample_id)
    assert_duplicate_metrics(actual)

    if actual != expected:
        expected_lines = json.dumps(expected, indent=2, sort_keys=True).splitlines()
        actual_lines = json.dumps(actual, indent=2, sort_keys=True).splitlines()
        print(
            "\n".join(
                difflib.unified_diff(
                    expected_lines,
                    actual_lines,
                    fromfile=str(expected_path),
                    tofile=str(args.root),
                    n=3,
                )
            )[:100_000]
        )
        raise AssertionError(f"Isolated Phase 5 {scenario} semantic regression mismatch")
    print(f"PASS: isolated {scenario} DNA family matches {expected_path}")


def validate_empty_coverage(args: argparse.Namespace) -> None:
    warning = (
        args.root
        / "pipeline_info/warnings/phase5_empty_Normal_H3K27ac.zero_mapped_nodup_bam.tsv"
    )
    if not warning.is_file():
        raise AssertionError(f"Empty-BAM warning was not published: {warning}")
    rows = list(csv.DictReader(warning.open(encoding="utf-8"), delimiter="\t"))
    if len(rows) != 1 or rows[0].get("mapped_reads") != "0":
        raise AssertionError(f"Unexpected empty-BAM warning: {rows}")
    if list((args.root / "dna_align").glob("phase5_empty*.bw")):
        raise AssertionError("Empty NoDup BAM unexpectedly produced coverage")


def main(args: argparse.Namespace) -> None:
    validate_trace(args)
    compare_scenario(args, "dna_single")
    compare_scenario(args, "dna_dual")
    validate_empty_coverage(args)
    print(f"PASS: Phase 5 {args.engine} execution preserved DNA family semantics")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--engine", required=True, choices=("docker", "conda"))
    parser.add_argument(
        "--expected-dir",
        type=Path,
        default=REPO / "tests/regression/contracts/v1.1.1",
    )
    parser.add_argument("--samtools", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
