#!/usr/bin/env python3

import argparse
import gzip
import os
import shutil
import subprocess
import sys
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


def extract_cr_and_others(comment: str):
    cb = ""
    um = ""
    others = []

    for token in comment.replace("\t", " ").split():
        if token.startswith("CB:"):
            cb = token.rsplit(":", 1)[-1]
        elif token.startswith("UM:"):
            um = token.rsplit(":", 1)[-1]
        else:
            others.append(token)

    if not cb or not um:
        return "", others

    return cb + um, others


def normalize_qname(name1: str, name2: str):
    if name1 == name2:
        return name1
    if len(name1) > 2 and len(name2) > 2 and name1.endswith("/1") and name2.endswith("/2") and name1[:-2] == name2[:-2]:
        return name1[:-2]
    return ""


def resolve_codon_bin() -> str:
    configured = os.environ.get("CODON_BIN")
    if configured:
        codon_bin = Path(configured)
        if not codon_bin.exists() or not os.access(codon_bin, os.X_OK):
            raise RuntimeError(f"Configured CODON_BIN is missing or not executable: {codon_bin}")
        return str(codon_bin)

    resolved = shutil.which("codon")
    if resolved is None:
        raise RuntimeError("codon executable not found in PATH")
    return resolved


def mock_fq_to_sam(args):
    with open(args.output_sam, "wt", encoding="utf-8") as out:
        out.write("@HD\tVN:1.6\tSO:unsorted\n")

        for r1_rec, r2_rec in zip(fastq_iter(args.r1), fastq_iter(args.r2)):
            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, _ = parse_header(r2_rec[0])

            if "NoMatch" in r1_comment:
                continue

            cr_value, other_tags = extract_cr_and_others(r1_comment)
            if not cr_value:
                raise ValueError(f"Missing CB or UM tag in FASTQ comment for {r1_name}")

            qname = normalize_qname(r1_name, r2_name)
            if not qname:
                raise ValueError(f"Mate names do not match: {r1_name} vs {r2_name}")

            out.write(
                f"{qname}\t77\t*\t0\t0\t*\t*\t0\t0\t{r1_rec[1]}\t{r1_rec[3]}\tCR:Z:{cr_value}"
            )
            for tag in other_tags:
                out.write(f"\t{tag}")
            out.write("\n")

            out.write(
                f"{qname}\t141\t*\t0\t0\t*\t*\t0\t0\t{r2_rec[1]}\t{r2_rec[3]}\tCR:Z:{cr_value}"
            )
            for tag in other_tags:
                out.write(f"\t{tag}")
            out.write("\n")


def real_fq_to_sam(args):
    codon_bin = resolve_codon_bin()

    cmd = [
        codon_bin,
        "run",
        "-plugin",
        "seq",
        "-release",
        str(args.script),
        str(args.r1),
        str(args.r2),
        str(args.output_sam),
    ]
    subprocess.run(cmd, check=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["real", "mock"])
    parser.add_argument("--script", required=True, type=Path)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--output-sam", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.mode == "mock":
        mock_fq_to_sam(args)
    else:
        real_fq_to_sam(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)
