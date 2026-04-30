#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from tresflow_fastq_utils import (
    fastq_iter,
    find_existing_output,
    log_event,
    open_maybe_gzip,
    parse_header,
    resolve_codon_bin,
    resolve_temp_root,
    strict_move_fastq,
    tagged_fastq_candidates,
)


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def mock_tag_umi(args):
    umi_counts = Counter()

    with open_maybe_gzip(args.output_r1, "wt") as out_r1, open_maybe_gzip(
        args.output_r2, "wt"
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
    codon_bin = resolve_codon_bin()

    with tempfile.TemporaryDirectory(prefix="tresflow_tag_umi_", dir=resolve_temp_root()) as tmpdir:
        tmp_path = Path(tmpdir)
        cmd = [
            codon_bin,
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
        codon_start = time.monotonic()
        log_event("Starting Codon Tag_UMI.codon", args.i2, args.r1, args.r2)
        subprocess.run(cmd, check=True)
        log_event("Finished Codon Tag_UMI.codon", args.i2, args.r1, args.r2, elapsed=time.monotonic() - codon_start)

        expected_r1 = find_existing_output(
            tmp_path,
            tagged_fastq_candidates(args.r1.name, args.tag),
            "tagged R1 FASTQ",
        )
        expected_r2 = find_existing_output(
            tmp_path,
            tagged_fastq_candidates(args.r2.name, args.tag),
            "tagged R2 FASTQ",
        )
        expected_counts = tmp_path / f"Reads_Per_Barcode_{args.sample}_{args.tag}.tsv"

        strict_move_fastq(expected_r1, args.output_r1)
        strict_move_fastq(expected_r2, args.output_r2)
        shutil.move(expected_counts, args.output_counts)

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
