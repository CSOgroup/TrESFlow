#!/usr/bin/env python3

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def copy_as_gzip(source: Path, destination: Path):
    with open_maybe_gzip(source, "rt") as src, gzip.open(destination, "wt") as dst:
        shutil.copyfileobj(src, dst)


def stem_without_fastq_suffix(name: str) -> str:
    if name.endswith(".fastq.gz"):
        return name[: -len(".fastq.gz")]
    if name.endswith(".fq.gz"):
        return name[: -len(".fq.gz")]
    if name.endswith(".fastq"):
        return name[: -len(".fastq")]
    if name.endswith(".fq"):
        return name[: -len(".fq")]
    return Path(name).stem


def mock_trim(args):
    copy_as_gzip(args.r1, args.output_r1)
    copy_as_gzip(args.r2, args.output_r2)


def real_trim(args):
    if shutil.which("trim_galore") is None:
        raise RuntimeError("trim_galore executable not found in PATH")

    with tempfile.TemporaryDirectory(prefix="tresflow_trim_galore_") as tmpdir:
        tmp_path = Path(tmpdir)
        cmd = [
            "trim_galore",
            "--quality",
            str(args.quality),
            "--cores",
            str(args.cores),
            "--output_dir",
            str(tmp_path),
            "--gzip",
            "--length",
            str(args.length),
            "--paired",
            str(args.r1),
            str(args.r2),
        ]
        subprocess.run(cmd, check=True)

        expected_r1 = tmp_path / f"{stem_without_fastq_suffix(args.r1.name)}_val_1.fq.gz"
        expected_r2 = tmp_path / f"{stem_without_fastq_suffix(args.r2.name)}_val_2.fq.gz"

        shutil.move(expected_r1, args.output_r1)
        shutil.move(expected_r2, args.output_r2)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["real", "mock"])
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--quality", required=True, type=int)
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--length", required=True, type=int)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.mode == "mock":
        mock_trim(args)
    else:
        real_trim(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)
