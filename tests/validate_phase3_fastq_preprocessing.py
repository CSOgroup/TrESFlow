#!/usr/bin/env python3
"""Semantic assertions for isolated Phase 3 Docker and Conda executions."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


EXPECTED_TASKS = Counter(
    {
        "TRIM_RNA_FASTQS": 1,
        "TRIM_DNA_FASTQS": 2,
        "DUAL_TAG_ARTIFACT_FILTER": 1,
        "COMPRESS_RNA_SPLIT_FASTQS": 1,
        "COMPRESS_DNA_SPLIT_FASTQS": 1,
    }
)
TRIM_CONTAINER = (
    "quay.io/biocontainers/trim-galore@"
    "sha256:a02bb87b8ce02d86efd0ffd65e2cce1559b52689faab42faad1df145657390cf"
)
CUTADAPT_CONTAINER = (
    "quay.io/biocontainers/cutadapt@"
    "sha256:2049f305574854edb189ccad7038fda4801ef16458bcde0239383d42d4a3f83a"
)
FORBIDDEN_TASK_COUPLING = (
    "/home/annan/micromamba/envs/tres",
    "${projectDir}",
    "PYTHON3_BIN",
    "TRIM_GALORE_BIN",
    "CUTADAPT_BIN",
    "PIGZ_BIN",
    "runtime.env_prefix",
)


def read_fastq(path: Path) -> list[tuple[str, str, str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    if len(lines) % 4:
        raise AssertionError(f"Incomplete FASTQ records: {path}")
    records = []
    for offset in range(0, len(lines), 4):
        record = tuple(lines[offset : offset + 4])
        if not record[0].startswith("@") or record[2] != "+":
            raise AssertionError(f"Malformed FASTQ record in {path}: {record!r}")
        if len(record[1]) != len(record[3]):
            raise AssertionError(f"Sequence/quality length mismatch in {path}: {record[0]}")
        records.append(record)
    return records


def read_trace(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def process_name(row: dict[str, str]) -> str:
    return row["name"].split(" (", 1)[0]


def sample_name(row: dict[str, str]) -> str:
    return row["name"].split(" (", 1)[1].rstrip(")")


def require_uncompressed_fastq(path: Path, expected_pairs: int) -> None:
    if path.read_bytes().startswith(b"\x1f\x8b"):
        raise AssertionError(f"Expected uncompressed FASTQ: {path}")
    if len(read_fastq(path)) != expected_pairs:
        raise AssertionError(f"Unexpected pair count in {path}")


def validate(args: argparse.Namespace) -> None:
    rows = read_trace(args.trace)
    observed = Counter(process_name(row) for row in rows)
    if observed != EXPECTED_TASKS:
        raise AssertionError(f"Unexpected Phase 3 task set: {observed}")
    for row in rows:
        if row["status"] != "COMPLETED" or row["exit"] != "0":
            raise AssertionError(f"Task did not complete: {row}")

        workdir = Path(row["workdir"])
        script = (workdir / ".command.sh").read_text(encoding="utf-8")
        wrapper = (workdir / ".command.run").read_text(encoding="utf-8")
        for token in FORBIDDEN_TASK_COUPLING:
            if token in script:
                raise AssertionError(f"Forbidden host coupling {token!r}: {row['name']}")
        if process_name(row).startswith("TRIM_"):
            if 'python3 "tresflow/bin/run_trim_galore.py"' not in script:
                raise AssertionError(f"Trim helper is not staged: {row['name']}")
        elif process_name(row) == "DUAL_TAG_ARTIFACT_FILTER":
            if 'python3 "tresflow/bin/run_dual_tag_artifact_filter.py"' not in script:
                raise AssertionError("Dual-tag helper is not staged")
        else:
            if "pigz -c -p" not in script:
                raise AssertionError(f"Compression does not use PATH pigz: {row['name']}")

        expected_container = (
            TRIM_CONTAINER if process_name(row).startswith("TRIM_") else CUTADAPT_CONTAINER
        )
        if row["container"] != expected_container:
            raise AssertionError(
                f"Unexpected declared container for {row['name']}: {row['container']}"
            )
        if args.engine == "docker":
            if "docker run " not in wrapper or expected_container not in wrapper:
                raise AssertionError(f"Task was not executed by Docker: {row['name']}")
        else:
            if "# conda environment" not in wrapper or "micromamba activate " not in wrapper:
                raise AssertionError(f"Task did not activate its Conda environment: {row['name']}")
            for executor in ("docker run ", "apptainer exec ", "singularity exec "):
                if executor in wrapper:
                    raise AssertionError(
                        f"Conda task unexpectedly used {executor.strip()}: {row['name']}"
                    )

    for row in rows:
        workdir = Path(row["workdir"])
        name = process_name(row)
        if name == "TRIM_RNA_FASTQS":
            require_uncompressed_fastq(
                workdir / "phase3_rna.sample_barcode_umi_cell.R1_val_1.fq", 3
            )
            require_uncompressed_fastq(
                workdir / "phase3_rna.sample_barcode_umi_cell.R2_val_2.fq", 3
            )
        elif name == "TRIM_DNA_FASTQS":
            sample = sample_name(row)
            require_uncompressed_fastq(
                workdir / f"{sample}.dna_sample_barcode_modality_cell.R1_val_1.fq", 3
            )
            require_uncompressed_fastq(
                workdir / f"{sample}.dna_sample_barcode_modality_cell.R2_val_2.fq", 3
            )
        elif name == "DUAL_TAG_ARTIFACT_FILTER":
            r1 = workdir / "phase3_dna_dual.dna_dual_tag_clean.R1.fastq"
            r2 = workdir / "phase3_dna_dual.dna_dual_tag_clean.R2.fastq"
            if [record[0] for record in read_fastq(r1)] != ["@salvage/1", "@clean/1"]:
                raise AssertionError("Dual-tag R1 filtering behavior changed")
            if [record[0] for record in read_fastq(r2)] != ["@salvage/2", "@clean/2"]:
                raise AssertionError("Dual-tag R2 filtering behavior changed")
            summary_path = workdir / "phase3_dna_dual.dual_tag_artifact_filter.summary.tsv"
            summary_rows = list(csv.DictReader(summary_path.open(encoding="utf-8"), delimiter="\t"))
            expected = {"input_pairs": "3", "retained_pairs": "2", "rejected_pairs": "1"}
            if len(summary_rows) != 1 or any(summary_rows[0][key] != value for key, value in expected.items()):
                raise AssertionError(f"Unexpected dual-tag summary: {summary_rows}")

    source_by_mate = {"R1": args.fixture_r1.read_bytes(), "R2": args.fixture_r2.read_bytes()}
    for modality in ("rna", "dna"):
        for mate in ("R1", "R2"):
            compressed = args.outdir / f"{modality}_split_fastqs" / f"phase3_Normal_{mate}.fastq.gz"
            with gzip.open(compressed, "rb") as handle:
                observed_bytes = handle.read()
            if observed_bytes != source_by_mate[mate]:
                raise AssertionError(f"Compressed publication changed FASTQ bytes: {compressed}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", required=True, choices=("docker", "conda"))
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--fixture-r1", required=True, type=Path)
    parser.add_argument("--fixture-r2", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    validate(parse_args())
