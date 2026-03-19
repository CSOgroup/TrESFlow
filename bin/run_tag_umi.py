#!/usr/bin/env python3

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def fastq_iter(path: Path):
    with open_maybe_gzip(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def parse_header(header: str):
    if not header.startswith("@"):
        raise ValueError(f"FASTQ header does not start with '@': {header}")
    body = header[1:]
    parts = body.split(" ", 1)
    if len(parts) == 1:
        return parts[0], ""
    return parts[0], parts[1]


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def mock_tag_umi(args):
    umi_counts = Counter()

    with open(args.output_r1, "wt", encoding="utf-8") as out_r1, open(
        args.output_r2, "wt", encoding="utf-8"
    ) as out_r2:
        for i2_rec, r1_rec, r2_rec in zip(
            fastq_iter(args.i2),
            fastq_iter(args.r1),
            fastq_iter(args.r2),
        ):
            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, _ = parse_header(r2_rec[0])
            if r1_name != r2_name:
                raise ValueError(f"Read name mismatch: {r1_name} != {r2_name}")

            umi = i2_rec[1][args.bc_start : args.bc_start + args.bc_len]
            if len(umi) != args.bc_len:
                raise ValueError(f"UMI slice shorter than expected in {args.i2}: {umi}")
            if args.rev_comp:
                umi = revcomp(umi)
            umi_counts[umi] += 1

            comment = f"{args.tag}:Z:{umi}"
            if r1_comment:
                comment = f"{comment}\t{r1_comment}"
            header = f"@{r1_name} {comment}"

            out_r1.write(f"{header}\n{r1_rec[1]}\n+\n{r1_rec[3]}\n")
            out_r2.write(f"{header}\n{r2_rec[1]}\n+\n{r2_rec[3]}\n")

    with open(args.output_counts, "wt", encoding="utf-8") as handle:
        for umi, count in umi_counts.items():
            handle.write(f"{count}\t{umi}\n")


def real_tag_umi(args):
    if shutil.which("codon") is None:
        raise RuntimeError("codon executable not found in PATH")

    with tempfile.TemporaryDirectory(prefix="tresflow_tag_umi_") as tmpdir:
        tmp_path = Path(tmpdir)
        cmd = [
            "codon",
            "run",
            "-plugin",
            "seq",
            "-release",
            "-D",
            f"BC_LEN={args.bc_len}",
            "-D",
            f"BC_START={args.bc_start}",
            str(args.script),
            str(args.i2),
            str(args.r1),
            str(args.r2),
            args.sample,
            args.tag,
            str(tmp_path),
        ]
        subprocess.run(cmd, check=True)

        expected_r1 = tmp_path / f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fastq"
        expected_r2 = tmp_path / f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fastq"
        expected_counts = tmp_path / f"Reads_Per_Barcode_{args.sample}_{args.tag}.tsv"

        shutil.move(expected_r1, args.output_r1)
        shutil.move(expected_r2, args.output_r2)
        shutil.move(expected_counts, args.output_counts)


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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["real", "mock"])
    parser.add_argument("--script", required=True, type=Path)
    parser.add_argument("--i2", required=True, type=Path)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--bc-len", required=True, type=int)
    parser.add_argument("--bc-start", required=True, type=int)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--output-counts", required=True, type=Path)
    parser.add_argument("--rev-comp", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.mode == "mock":
        mock_tag_umi(args)
    else:
        real_tag_umi(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)
