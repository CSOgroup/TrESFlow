#!/usr/bin/env python3
"""Discard dual-tagmentation read pairs containing audited linker signatures."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


EXPECTED_SIGNATURE_SHA256 = "67c6f1789ef5e36492562203ac38fc13fa901058047ed2bd37b304d85a30ae0f"
EXPECTED_SIGNATURE_COUNT = 48
SIGNATURE_LENGTH = 23
SUMMARY_COLUMNS = (
    "sample_id",
    "tagmentation",
    "input_pairs",
    "retained_pairs",
    "rejected_pairs",
    "rejected_fraction",
    "r1_with_signature",
    "r2_with_signature",
    "cutadapt_version",
    "signature_fasta_sha256",
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_and_validate_signatures(path: Path) -> list[str]:
    observed_hash = file_sha256(path)
    if observed_hash != EXPECTED_SIGNATURE_SHA256:
        raise ValueError(
            f"Unexpected dual-tag signature FASTA SHA-256 for {path}: "
            f"expected {EXPECTED_SIGNATURE_SHA256}, observed {observed_hash}"
        )

    sequences = []
    current = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                sequences.append("".join(current).upper())
                current = []
        else:
            current.append(line)
    if current:
        sequences.append("".join(current).upper())

    if len(sequences) != EXPECTED_SIGNATURE_COUNT:
        raise ValueError(
            f"Dual-tag signature FASTA must contain {EXPECTED_SIGNATURE_COUNT} sequences; "
            f"observed {len(sequences)}"
        )
    if len(set(sequences)) != EXPECTED_SIGNATURE_COUNT:
        raise ValueError("Dual-tag signature FASTA contains duplicate sequences")
    for sequence in sequences:
        if len(sequence) != SIGNATURE_LENGTH or any(base not in "ACGT" for base in sequence):
            raise ValueError(
                f"Dual-tag signatures must contain exactly {SIGNATURE_LENGTH} A/C/G/T bases: {sequence}"
            )
    return sequences


def resolve_cutadapt_bin(configured: Path | None) -> Path:
    if configured is None:
        resolved = shutil.which("cutadapt")
        if resolved is None:
            raise RuntimeError("cutadapt executable not found in PATH")
        configured = Path(resolved)
    if not configured.exists() or not os.access(configured, os.X_OK):
        raise RuntimeError(f"cutadapt executable not found or not executable: {configured}")
    return configured


def build_cutadapt_command(args) -> list[str]:
    cutadapt_bin = resolve_cutadapt_bin(args.cutadapt_bin)
    return [
        str(cutadapt_bin),
        "-j",
        str(args.cpus),
        "-e",
        "0",
        "-O",
        str(SIGNATURE_LENGTH),
        "--no-indels",
        "--action=none",
        "-b",
        f"file:{args.signature_fasta}",
        "-B",
        f"file:{args.signature_fasta}",
        "--discard-trimmed",
        "--pair-filter=any",
        "--json",
        str(args.output_json),
        "-o",
        str(args.output_r1),
        "-p",
        str(args.output_r2),
        str(args.r1),
        str(args.r2),
    ]


def _required_nonnegative_integer(mapping: dict, key: str) -> int:
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"Cutadapt JSON field '{key}' must be a non-negative integer: {value!r}")
    return value


def summary_from_cutadapt_json(
    sample_id: str,
    tagmentation: str,
    json_path: Path,
    signature_fasta: Path,
) -> dict:
    report = json.loads(json_path.read_text(encoding="utf-8"))
    read_counts = report.get("read_counts")
    if not isinstance(read_counts, dict):
        raise ValueError("Cutadapt JSON is missing the read_counts object")
    filtered = read_counts.get("filtered")
    if not isinstance(filtered, dict):
        raise ValueError("Cutadapt JSON is missing read_counts.filtered")

    input_pairs = _required_nonnegative_integer(read_counts, "input")
    retained_pairs = _required_nonnegative_integer(read_counts, "output")
    rejected_pairs = _required_nonnegative_integer(filtered, "discard_trimmed")
    r1_with_signature = _required_nonnegative_integer(read_counts, "read1_with_adapter")
    r2_with_signature = _required_nonnegative_integer(read_counts, "read2_with_adapter")
    if input_pairs != retained_pairs + rejected_pairs:
        raise ValueError(
            "Cutadapt pair accounting mismatch: "
            f"input={input_pairs}, retained={retained_pairs}, rejected={rejected_pairs}"
        )
    cutadapt_version = str(report.get("cutadapt_version") or "").strip()
    if not cutadapt_version:
        raise ValueError("Cutadapt JSON is missing cutadapt_version")

    return {
        "sample_id": sample_id,
        "tagmentation": tagmentation,
        "input_pairs": input_pairs,
        "retained_pairs": retained_pairs,
        "rejected_pairs": rejected_pairs,
        "rejected_fraction": (rejected_pairs / input_pairs) if input_pairs else 0.0,
        "r1_with_signature": r1_with_signature,
        "r2_with_signature": r2_with_signature,
        "cutadapt_version": cutadapt_version,
        "signature_fasta_sha256": file_sha256(signature_fasta),
    }


def write_summary(path: Path, summary: dict) -> None:
    values = []
    for column in SUMMARY_COLUMNS:
        value = summary[column]
        if column == "rejected_fraction":
            value = f"{float(value):.8f}"
        values.append(str(value))
    path.write_text("\t".join(SUMMARY_COLUMNS) + "\n" + "\t".join(values) + "\n", encoding="utf-8")


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")


def fastq_records(path: Path):
    with open_maybe_gzip(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not quality:
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header, sequence, plus, quality


def mock_filter(args, signatures: list[str]) -> None:
    counts = {"input": 0, "output": 0, "discard_trimmed": 0, "r1": 0, "r2": 0}
    with args.output_r1.open("wt", encoding="utf-8") as output_r1, args.output_r2.open(
        "wt", encoding="utf-8"
    ) as output_r2:
        r1_iterator = fastq_records(args.r1)
        r2_iterator = fastq_records(args.r2)
        while True:
            r1_record = next(r1_iterator, None)
            r2_record = next(r2_iterator, None)
            if r1_record is None and r2_record is None:
                break
            if r1_record is None or r2_record is None:
                raise ValueError("Paired FASTQ inputs contain different numbers of records")
            counts["input"] += 1
            r1_match = any(signature in r1_record[1].strip().upper() for signature in signatures)
            r2_match = any(signature in r2_record[1].strip().upper() for signature in signatures)
            counts["r1"] += int(r1_match)
            counts["r2"] += int(r2_match)
            if r1_match or r2_match:
                counts["discard_trimmed"] += 1
                continue
            counts["output"] += 1
            output_r1.writelines(r1_record)
            output_r2.writelines(r2_record)

    report = {
        "tag": "Cutadapt report",
        "schema_version": [0, 3],
        "cutadapt_version": "mock",
        "read_counts": {
            "input": counts["input"],
            "filtered": {"discard_trimmed": counts["discard_trimmed"]},
            "output": counts["output"],
            "read1_with_adapter": counts["r1"],
            "read2_with_adapter": counts["r2"],
        },
    }
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def run(args) -> None:
    signatures = load_and_validate_signatures(args.signature_fasta)
    if args.cpus < 1:
        raise ValueError("--cpus must be at least 1")
    if args.tagmentation != "dual":
        raise ValueError("The dual-tag artifact filter may only run with --tagmentation dual")

    if args.mode == "mock":
        mock_filter(args, signatures)
    else:
        subprocess.run(build_cutadapt_command(args), check=True)

    summary = summary_from_cutadapt_json(
        args.sample_id,
        args.tagmentation,
        args.output_json,
        args.signature_fasta,
    )
    write_summary(args.output_summary, summary)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=("real", "mock"))
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--tagmentation", required=True, choices=("dual",))
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--signature-fasta", required=True, type=Path)
    parser.add_argument("--cutadapt-bin", type=Path)
    parser.add_argument("--cpus", required=True, type=int)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def main():
    run(parse_args())


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        sys.stderr.write(f"ERROR: {error}\n")
        sys.exit(1)
