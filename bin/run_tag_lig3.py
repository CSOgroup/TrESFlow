#!/usr/bin/env python3

import argparse
import gzip
import os
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


def open_maybe_gzip_binary(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def move_or_compress_fastq(source: Path, destination: Path):
    if source.suffix == destination.suffix:
        shutil.move(source, destination)
        return

    with open_maybe_gzip_binary(source, "rb") as src, open_maybe_gzip_binary(destination, "wb") as dst:
        shutil.copyfileobj(src, dst)
    source.unlink()


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


def load_whitelist(path: Path):
    with open(path, "rt", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


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


def find_tag_value(comment: str, tag_name: str):
    for token in comment.replace("\t", " ").split():
        if token.startswith(f"{tag_name}:"):
            return token.rsplit(":", 1)[-1]
    return ""


def write_stats(output_stats: Path, n_reads: int, bc_reads: int, mismatch_stats, hd: int):
    reads_without_bc = n_reads - bc_reads
    n_passing_reads = bc_reads
    reads_in_too_small_bc = bc_reads - n_passing_reads

    lines = [
        f"reads\t{n_reads}\t{percent(n_reads, n_reads)}",
        f"reads_with_ligation\t{bc_reads}\t{percent(bc_reads, n_reads)}",
        f"reads_without_ligation\t{reads_without_bc}\t{percent(reads_without_bc, n_reads)}",
    ]
    for idx in range(hd + 1):
        lines.append(
            f"reads_with_ligation_at_{idx}_mismatches\t{mismatch_stats[idx]}\t{percent(mismatch_stats[idx], n_reads)}"
        )
    lines.extend(
        [
            f"reads_with_ligation_all_segments\t{n_passing_reads}\t{percent(n_passing_reads, n_reads)}",
            f"reads_with_ligation_missing_other_ligation\t{reads_in_too_small_bc}\t{percent(reads_in_too_small_bc, n_reads)}",
        ]
    )
    output_stats.write_text("\n".join(lines) + "\n", encoding="utf-8")


def percent(count: int, total: int) -> str:
    if total == 0:
        return "0.0%"
    return f"{(count / total) * 100.0}%"


def mock_tag(args):
    whitelist = load_whitelist(args.whitelist)
    start_positions = [15, 53, 91]
    lig_stats = [[0 for _ in range(args.hd + 1)] for _ in range(3)]
    barcode_counts = Counter()
    total_reads = 0

    with open_maybe_gzip(args.output_r1, "wt") as out_r1, \
        open_maybe_gzip(args.output_r2, "wt") as out_r2, \
        open(args.output_tag_records, "wt", encoding="utf-8") as tag_records:
        for i1_rec, r1_rec, r2_rec in zip(
            fastq_iter(args.i1),
            fastq_iter(args.r1),
            fastq_iter(args.r2),
        ):
            total_reads += 1

            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, _ = parse_header(r2_rec[0])
            i1_name, _ = parse_header(i1_rec[0])
            if len({i1_name, r1_name, r2_name}) != 1:
                raise ValueError(f"Read name mismatch: {i1_name}, {r1_name}, {r2_name}")

            sb_tag = find_tag_value(r1_comment, "SB")
            if not sb_tag:
                raise ValueError(f"Missing SB tag in FASTQ comment for {r1_name}")

            ligation_values = []
            for idx, start in enumerate(start_positions):
                barcode = i1_rec[1][start : start + args.bc_len]
                if len(barcode) != args.bc_len:
                    raise ValueError(f"Ligation barcode slice shorter than expected in {args.i1}: {barcode}")
                if barcode in whitelist:
                    lig_stats[idx][0] += 1
                    ligation_values.append(barcode)
                else:
                    ligation_values.append("NoMatch")

            final_bc = "NoMatch" if "NoMatch" in ligation_values else sb_tag + "".join(ligation_values)
            barcode_counts[final_bc] += 1

            header_comment = f"{args.tag}:Z:{final_bc}\tRG:Z:{final_bc}\t{r1_comment}"
            header = f"@{r1_name} {header_comment}"

            out_r1.write(f"{header}\n{r1_rec[1]}\n+\n{r1_rec[3]}\n")
            out_r2.write(f"{header}\n{r2_rec[1]}\n+\n{r2_rec[3]}\n")

            tag_records.write(
                f"{r1_name}\t{args.tag}:Z:{final_bc}\tL1:Z:{ligation_values[0]}\tL2:Z:{ligation_values[1]}\tL3:Z:{ligation_values[2]}\t{r1_comment}\n"
            )

    with open(args.output_counts, "wt", encoding="utf-8") as handle:
        for barcode, count in barcode_counts.items():
            handle.write(f"{count}\t{barcode}\n")

    for idx in range(3):
        bc_reads = lig_stats[idx][0]
        write_stats(args.output_stats[idx], total_reads, bc_reads, lig_stats[idx], args.hd)


def real_tag(args):
    codon_bin = resolve_codon_bin()

    with tempfile.TemporaryDirectory(prefix="tresflow_tag_lig3_") as tmpdir:
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
            f"HD={args.hd}",
            str(args.script),
            str(args.i1),
            str(args.r1),
            str(args.r2),
            str(args.whitelist),
            args.sample,
            args.tag,
            str(tmp_path),
        ]
        subprocess.run(cmd, check=True)

        expected_r1 = find_existing_output(
            tmp_path,
            [
                f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fastq.gz",
                f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fq.gz",
                f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fastq",
                f"{stem_without_fastq_suffix(args.r1.name)}_{args.tag}.fq",
            ],
            "tagged R1 FASTQ",
        )
        expected_r2 = find_existing_output(
            tmp_path,
            [
                f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fastq.gz",
                f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fq.gz",
                f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fastq",
                f"{stem_without_fastq_suffix(args.r2.name)}_{args.tag}.fq",
            ],
            "tagged R2 FASTQ",
        )
        expected_counts = tmp_path / f"Reads_Per_Barcode_{args.sample}_{args.tag}.tsv"
        expected_tag_records = tmp_path / f"Tag_Records_{args.sample}.tsv"
        expected_stats = [
            tmp_path / f"Barcode_Statistics_{args.sample}_{args.tag}_L1.tsv",
            tmp_path / f"Barcode_Statistics_{args.sample}_{args.tag}_L2.tsv",
            tmp_path / f"Barcode_Statistics_{args.sample}_{args.tag}_L3.tsv",
        ]

        move_or_compress_fastq(expected_r1, args.output_r1)
        move_or_compress_fastq(expected_r2, args.output_r2)
        shutil.move(expected_counts, args.output_counts)
        shutil.move(expected_tag_records, args.output_tag_records)
        for expected, output in zip(expected_stats, args.output_stats):
            shutil.move(expected, output)


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
    parser.add_argument("--i1", required=True, type=Path)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--whitelist", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--bc-len", required=True, type=int)
    parser.add_argument("--hd", required=True, type=int)
    parser.add_argument("--output-r1", required=True, type=Path)
    parser.add_argument("--output-r2", required=True, type=Path)
    parser.add_argument("--output-counts", required=True, type=Path)
    parser.add_argument("--output-tag-records", required=True, type=Path)
    parser.add_argument("--output-stats", action="append", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    if len(args.output_stats) != 3:
        raise ValueError("Expected exactly three --output-stats arguments")
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
