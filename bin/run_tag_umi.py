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
    fastq_input_spec,
    find_existing_output,
    log_event,
    open_maybe_gzip,
    parse_header,
    read_read_set_counts,
    resolve_codon_bin,
    resolve_fastq_paths,
    resolve_temp_root,
    strict_move_fastq,
    synchronized_fastq_iter,
    tagged_fastq_candidates,
    write_read_set_counts,
)


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def mock_tag_umi(args):
    umi_counts = Counter()
    total_reads = 0
    observed_read_set_counts = []

    with open_maybe_gzip(args.output_r1, "wt") as out_r1, open_maybe_gzip(
        args.output_r2, "wt"
    ) as out_r2:
        for record_number, (i2_rec, r1_rec, r2_rec) in enumerate(
            synchronized_fastq_iter(
                {"i2": args.i2_paths, "r1": args.r1_paths, "r2": args.r2_paths},
                args.read_set_counts_values,
                observed_read_set_counts,
            ),
            start=1,
        ):
            r1_name, r1_comment = parse_header(r1_rec[0])

            umi = i2_rec[1][args.bc_start : args.bc_start + args.bc_len]
            if len(umi) != args.bc_len:
                raise ValueError(
                    f"UMI slice shorter than expected in virtual i2 stream at read {record_number}: {umi}"
                )
            if args.rev_comp:
                umi = revcomp(umi)
            umi_counts[umi] += 1

            comment = f"{args.tag}:Z:{umi}"
            if r1_comment:
                comment = f"{comment}\t{r1_comment}"
            header = f"@{r1_name} {comment}"

            out_r1.write(f"{header}\n{r1_rec[1]}\n+\n{r1_rec[3]}\n")
            out_r2.write(f"{header}\n{r2_rec[1]}\n+\n{r2_rec[3]}\n")
            total_reads += 1

    if total_reads == 0:
        raise ValueError("Synchronized FASTQ streams contain no records")
    if args.output_read_set_counts is not None:
        write_read_set_counts(args.output_read_set_counts, observed_read_set_counts)

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
            fastq_input_spec(args.i2, args.i2_manifest),
            fastq_input_spec(args.r1, args.r1_manifest),
            fastq_input_spec(args.r2, args.r2_manifest),
            args.sample,
            args.tag,
            str(tmp_path),
            "tresflow_tag_umi.R1.fastq",
            "tresflow_tag_umi.R2.fastq",
            str(args.read_set_counts.resolve()) if args.read_set_counts is not None else "-",
        ]
        codon_start = time.monotonic()
        log_event("Starting Codon Tag_UMI.codon", *(args.i2_paths + args.r1_paths + args.r2_paths))
        subprocess.run(cmd, check=True)
        log_event(
            "Finished Codon Tag_UMI.codon",
            *(args.i2_paths + args.r1_paths + args.r2_paths),
            elapsed=time.monotonic() - codon_start,
        )

        expected_r1 = find_existing_output(
            tmp_path,
            ["tresflow_tag_umi.R1.fastq"] + tagged_fastq_candidates(args.r1_paths[0].name, args.tag),
            "tagged R1 FASTQ",
        )
        expected_r2 = find_existing_output(
            tmp_path,
            ["tresflow_tag_umi.R2.fastq"] + tagged_fastq_candidates(args.r2_paths[0].name, args.tag),
            "tagged R2 FASTQ",
        )
        expected_counts = tmp_path / f"Reads_Per_Barcode_{args.sample}_{args.tag}.tsv"
        expected_read_set_counts = tmp_path / "technical_read_set_counts.tsv"

        strict_move_fastq(expected_r1, args.output_r1)
        strict_move_fastq(expected_r2, args.output_r2)
        shutil.move(expected_counts, args.output_counts)
        if args.output_read_set_counts is not None:
            shutil.move(expected_read_set_counts, args.output_read_set_counts)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["real", "mock"])
    parser.add_argument("--script", required=True, type=Path)
    i2_group = parser.add_mutually_exclusive_group(required=True)
    i2_group.add_argument("--i2", type=Path)
    i2_group.add_argument("--i2-manifest", type=Path)
    r1_group = parser.add_mutually_exclusive_group(required=True)
    r1_group.add_argument("--r1", type=Path)
    r1_group.add_argument("--r1-manifest", type=Path)
    r2_group = parser.add_mutually_exclusive_group(required=True)
    r2_group.add_argument("--r2", type=Path)
    r2_group.add_argument("--r2-manifest", type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--bc-len", required=True, type=int)
    parser.add_argument("--bc-start", required=True, type=int)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--output-counts", required=True, type=Path)
    parser.add_argument("--rev-comp", action="store_true")
    parser.add_argument("--read-set-counts", type=Path)
    parser.add_argument("--output-read-set-counts", type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    args.i2_paths = resolve_fastq_paths(args.i2, args.i2_manifest)
    args.r1_paths = resolve_fastq_paths(args.r1, args.r1_manifest)
    args.r2_paths = resolve_fastq_paths(args.r2, args.r2_manifest)
    args.read_set_counts_values = (
        read_read_set_counts(args.read_set_counts) if args.read_set_counts is not None else None
    )
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
