#!/usr/bin/env python3

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path


def resolve_temp_root() -> Path:
    configured = os.environ.get("TMPDIR")
    if not configured:
        raise RuntimeError("TMPDIR is not set. Configure runtime.tmpdir in the samplesheet.")
    root = Path(configured).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    return root.resolve()


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


def normalize_sb_drop_first(sb: str):
    if len(sb) < 2:
        raise ValueError(f"SB tag length < 2: {sb}")
    return sb[1:]


def find_tag_value(comment: str, tag_name: str):
    for token in comment.replace("\t", " ").split():
        if token.startswith(f"{tag_name}:"):
            return token.rsplit(":", 1)[-1]
    return ""


def load_sb_group_map(path: Path, sample: str):
    sb_to_group = {}
    group_names = []
    group_seen = set()

    with open(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            row_sample, group_name, sb_bc = parts[0], parts[1], parts[2]
            if row_sample != sample:
                continue

            if group_name not in group_seen:
                group_names.append(group_name)
                group_seen.add(group_name)

            if sb_bc in sb_to_group and sb_to_group[sb_bc] != group_name:
                raise ValueError(f"SB group conflict for sample {sample} SB {sb_bc}")
            sb_to_group[sb_bc] = group_name

    if not sb_to_group:
        raise ValueError(f"No SB group mapping found for sample {sample} in {path}")

    return sb_to_group, group_names


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


def resolve_group(sample: str, sb_raw: str, sb_to_group):
    if sb_raw in sb_to_group:
        return sb_to_group[sb_raw]

    key = normalize_sb_drop_first(sb_raw)
    if key in sb_to_group:
        return sb_to_group[key]

    raise ValueError(f"SB not found in SB group map for sample {sample}: raw={sb_raw} key={key}")


def write_fastq_record(handle, name: str, comment: str, seq: str, qual: str):
    handle.write("@")
    handle.write(name)
    if comment:
        handle.write(" ")
        handle.write(comment)
    handle.write("\n")
    handle.write(seq)
    handle.write("\n+\n")
    handle.write(qual)
    handle.write("\n")


def write_rg_header(path: Path, sample: str, library_name: str, barcodes):
    with open(path, "wt", encoding="utf-8") as handle:
        for barcode in sorted(barcodes):
            handle.write(f"@RG\tID:{barcode}\tSM:{sample}\tLB:{library_name}\tPL:ELEMENT\tPM:AVITI_500MIO\n")


def mock_split(args):
    sb_to_group, group_names = load_sb_group_map(args.sb_group_map, args.sample)

    r1_handles = OrderedDict()
    r2_handles = OrderedDict()
    group_barcodes = {}

    for group_name in group_names:
        r1_handles[group_name] = gzip.open(args.output_dir / f"{args.sample}_{group_name}_R1.fq.gz", "wt")
        r2_handles[group_name] = gzip.open(args.output_dir / f"{args.sample}_{group_name}_R2.fq.gz", "wt")
        group_barcodes[group_name] = set()

    try:
        for r1_rec, r2_rec in zip(fastq_iter(args.r1), fastq_iter(args.r2)):
            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, r2_comment = parse_header(r2_rec[0])
            if r1_name != r2_name:
                raise ValueError(f"Read name mismatch: {r1_name} != {r2_name}")

            if "NoMatch" in r1_comment:
                continue

            cb = find_tag_value(r1_comment, "CB")
            sb = find_tag_value(r1_comment, "SB")
            if not cb or not sb:
                raise ValueError(f"Missing CB or SB tag in FASTQ comment for {r1_name}")

            group_name = resolve_group(args.sample, sb, sb_to_group)
            write_fastq_record(r1_handles[group_name], r1_name, r1_comment, r1_rec[1], r1_rec[3])
            write_fastq_record(r2_handles[group_name], r2_name, r2_comment, r2_rec[1], r2_rec[3])
            group_barcodes[group_name].add(cb)
    finally:
        for handle in r1_handles.values():
            handle.close()
        for handle in r2_handles.values():
            handle.close()

    for group_name in group_names:
        header_path = args.output_dir / f"SAM_RG_Header_{args.sample}_{group_name}.tsv"
        write_rg_header(header_path, args.sample, args.library_name, group_barcodes[group_name])


def real_split(args):
    codon_bin = resolve_codon_bin()

    with tempfile.TemporaryDirectory(prefix="tresflow_split_reads_rna_", dir=resolve_temp_root()) as tmpdir:
        tmp_path = Path(tmpdir)
        cmd = [
            codon_bin,
            "run",
            "-plugin",
            "seq",
            "-release",
            str(args.script),
            args.sample,
            str(tmp_path),
            args.library_name,
            "rna",
            "-",
            str(args.r1),
            str(args.r2),
            str(args.sb_group_map),
        ]
        subprocess.run(cmd, check=True)

        moved = 0
        for pattern in (
            f"{args.sample}_*_R1.fq.gz",
            f"{args.sample}_*_R2.fq.gz",
            f"SAM_RG_Header_{args.sample}_*.tsv",
        ):
            for source in sorted(tmp_path.glob(pattern)):
                shutil.move(source, args.output_dir / source.name)
                moved += 1

        if moved == 0:
            raise RuntimeError(f"No RNA split outputs were produced for sample {args.sample}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["real", "mock"])
    parser.add_argument("--script", required=True, type=Path)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--sb-group-map", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--library-name", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.mode == "mock":
        mock_split(args)
    else:
        real_split(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)
