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


def load_whitelist(path: Path):
    with open(path, "rt", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def load_whitelist_from_sb_group_map(path: Path, sample: str):
    sb_to_group = {}

    with open(path, "rt", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                raise ValueError(
                    f"Malformed sample-barcode group map line {line_no} in {path}: expected at least 3 columns"
                )

            row_sample, group_name, sb_bc = parts[0], parts[1], parts[2]
            if (
                row_sample.lower() == "sample"
                and group_name.lower() == "sb_group"
                and sb_bc.lower() == "sb_bc"
            ):
                continue

            if row_sample != sample:
                continue

            if sb_bc in sb_to_group and sb_to_group[sb_bc] != group_name:
                raise ValueError(
                    f"Conflicting sample-barcode group map rows for sample {sample} barcode {sb_bc}: "
                    f"{sb_to_group[sb_bc]} vs {group_name}"
                )

            sb_to_group[sb_bc] = group_name

    if not sb_to_group:
        raise ValueError(f"No sample-barcode group map rows found in {path} for sample {sample}")

    return set(sb_to_group.keys())


def write_stats(output_stats: Path, n_reads: int, bc_reads: int, mismatch_stats, hd: int):
    reads_without_bc = n_reads - bc_reads
    lines = [
        f"reads\t{n_reads}\t{percent(n_reads, n_reads)}",
        f"bc_reads\t{bc_reads}\t{percent(bc_reads, n_reads)}",
        f"reads_without_bc\t{reads_without_bc}\t{percent(reads_without_bc, n_reads)}",
    ]
    for idx in range(hd + 1):
        lines.append(
            f"bc_reads_with_{idx}_mismatches\t{mismatch_stats[idx]}\t{percent(mismatch_stats[idx], n_reads)}"
        )
    output_stats.write_text("\n".join(lines) + "\n", encoding="utf-8")


def percent(count: int, total: int) -> str:
    if total == 0:
        return "0.0%"
    return f"{(count / total) * 100.0}%"


def resolve_whitelist_values(args):
    if args.sb_group_map is not None:
        return load_whitelist_from_sb_group_map(args.sb_group_map, args.sample)
    return load_whitelist(args.whitelist)


def mock_tag(args):
    whitelist = resolve_whitelist_values(args)
    mismatch_stats = [0 for _ in range(args.hd + 1)]
    barcode_counts = Counter()
    total_reads = 0

    if args.first_pass_arg == "first_pass":
        first_pass = True
        prefix = ""
    elif args.first_pass_arg.startswith("first_pass_withBC_"):
        first_pass = True
        prefix = args.first_pass_arg[-1]
    elif args.first_pass_arg == "not_first_pass":
        first_pass = False
        prefix = ""
    else:
        raise ValueError(f"Unsupported first_pass_arg for mock mode: {args.first_pass_arg}")

    if args.rev_comp_arg not in {"rev", "fw"}:
        raise ValueError(f"Unsupported rev_comp_arg for mock mode: {args.rev_comp_arg}")

    with open(args.output_r1, "wt", encoding="utf-8") as out_r1, open(
        args.output_r2, "wt", encoding="utf-8"
    ) as out_r2:
        for i2_rec, r1_rec, r2_rec in zip(
            fastq_iter(args.i2),
            fastq_iter(args.r1),
            fastq_iter(args.r2),
        ):
            total_reads += 1
            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, _ = parse_header(r2_rec[0])
            if r1_name != r2_name:
                raise ValueError(f"Read name mismatch: {r1_name} != {r2_name}")

            barcode = i2_rec[1][args.bc_start : args.bc_start + args.bc_len]
            if len(barcode) != args.bc_len:
                raise ValueError(f"Barcode slice shorter than expected in {args.i2}: {barcode}")
            if args.rev_comp_arg == "rev":
                barcode = revcomp(barcode)

            if barcode in whitelist:
                mismatch_stats[0] += 1
                barcode_counts[barcode] += 1
                tag_value = prefix + barcode
            else:
                tag_value = "NoMatch"
                barcode_counts["NoMatch"] += 1

            comment_parts = [f"{args.tag}:Z:{tag_value}"]
            if (not first_pass) and r1_comment:
                comment_parts.append(r1_comment)
            comment = "\t".join(comment_parts)
            header = f"@{r1_name} {comment}".rstrip()

            out_r1.write(f"{header}\n{r1_rec[1]}\n+\n{r1_rec[3]}\n")
            out_r2.write(f"{header}\n{r2_rec[1]}\n+\n{r2_rec[3]}\n")

    with open(args.output_counts, "wt", encoding="utf-8") as handle:
        for barcode, count in barcode_counts.items():
            handle.write(f"{count}\t{barcode}\n")

    write_stats(args.output_stats, total_reads, mismatch_stats[0], mismatch_stats, args.hd)


def real_tag(args):
    if shutil.which("codon") is None:
        raise RuntimeError("codon executable not found in PATH")

    with tempfile.TemporaryDirectory(prefix="tresflow_tag_") as tmpdir:
        tmp_path = Path(tmpdir)
        if args.sb_group_map is not None:
            whitelist_path = tmp_path / f"{args.sample}_{args.tag}.derived_whitelist.txt"
            whitelist_values = sorted(load_whitelist_from_sb_group_map(args.sb_group_map, args.sample))
            whitelist_path.write_text("\n".join(whitelist_values) + "\n", encoding="utf-8")
        else:
            whitelist_path = args.whitelist

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
            "-D",
            f"HD={args.hd}",
            str(args.script),
            str(args.i2),
            str(args.r1),
            str(args.r2),
            str(whitelist_path),
            args.sample,
            args.tag,
            str(tmp_path),
            args.first_pass_arg,
            args.rev_comp_arg,
        ]
        subprocess.run(cmd, check=True)

        expected_r1 = find_existing_output(
            tmp_path,
            [
                f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fastq",
                f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fq",
            ],
            "tagged R1 FASTQ",
        )
        expected_r2 = find_existing_output(
            tmp_path,
            [
                f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fastq",
                f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fq",
            ],
            "tagged R2 FASTQ",
        )
        expected_counts = tmp_path / f"Reads_Per_Barcode_{args.sample}_{args.tag}.tsv"
        expected_stats = tmp_path / f"Barcode_Statistics_{args.sample}_{args.tag}.tsv"

        shutil.move(expected_r1, args.output_r1)
        shutil.move(expected_r2, args.output_r2)
        shutil.move(expected_counts, args.output_counts)
        shutil.move(expected_stats, args.output_stats)


def find_existing_output(base_dir: Path, candidate_names, label: str) -> Path:
    for candidate_name in candidate_names:
        candidate = base_dir / candidate_name
        if candidate.exists():
            return candidate
    joined = ", ".join(str(base_dir / name) for name in candidate_names)
    raise FileNotFoundError(f"Expected {label} in one of: {joined}")


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
    whitelist_group = parser.add_mutually_exclusive_group(required=True)
    whitelist_group.add_argument("--sb-group-map", type=Path)
    whitelist_group.add_argument("--whitelist", type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--bc-len", required=True, type=int)
    parser.add_argument("--bc-start", required=True, type=int)
    parser.add_argument("--hd", required=True, type=int)
    parser.add_argument("--first-pass-arg", required=True)
    parser.add_argument("--rev-comp-arg", required=True)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--output-counts", required=True, type=Path)
    parser.add_argument("--output-stats", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.mode == "mock":
        mock_tag(args)
    else:
        real_tag(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)
