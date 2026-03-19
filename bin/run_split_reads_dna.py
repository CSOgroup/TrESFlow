#!/usr/bin/env python3

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
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


def load_mo_map(path: Path, sample: str, group_names):
    group_to_index = {group_name: idx for idx, group_name in enumerate(group_names)}
    mark_names = []
    mark_seen = set()
    mappings = []
    pair_to_mark = {}
    seen_format = 0

    with open(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            row_sample = parts[0]
            if row_sample != sample:
                continue

            if len(parts) >= 4:
                if seen_format == 0:
                    seen_format = 4
                elif seen_format != 4:
                    raise ValueError(f"Mixed 3-col and 4-col MO map lines for sample {sample}")

                group_name, mark_name, mo_bc = parts[1], parts[2], parts[3]
                if group_name not in group_to_index:
                    raise ValueError(
                        f"MO map uses sb_group '{group_name}' not present in SB group map for sample {sample}"
                    )
                group_index = group_to_index[group_name]
            else:
                if seen_format == 0:
                    seen_format = 3
                elif seen_format != 3:
                    raise ValueError(f"Mixed 3-col and 4-col MO map lines for sample {sample}")

                mark_name, mo_bc = parts[1], parts[2]
                group_index = -1

            if mark_name not in mark_seen:
                mark_names.append(mark_name)
                mark_seen.add(mark_name)

            key = (group_index, mo_bc)
            if key in pair_to_mark and pair_to_mark[key] != mark_name:
                raise ValueError(
                    f"MO map conflict for sample {sample} gid {group_index} MO {mo_bc}: "
                    f"{pair_to_mark[key]} vs {mark_name}"
                )

            if key not in pair_to_mark:
                pair_to_mark[key] = mark_name
                mappings.append((group_index, mo_bc, mark_name))

    if not mappings:
        raise ValueError(f"No MO mapping found for sample {sample} in {path}")

    return mappings, mark_names


def resolve_group(sample: str, sb_raw: str, sb_to_group):
    if sb_raw in sb_to_group:
        return sb_to_group[sb_raw]

    key = normalize_sb_drop_first(sb_raw)
    if key in sb_to_group:
        return sb_to_group[key]

    raise ValueError(f"SB not found in SB group map for sample {sample}: raw={sb_raw} key={key}")


def find_mark_for_mo(mo: str, group_name: str, group_names, mappings):
    group_index = group_names.index(group_name)

    for gid, mo_bc, mark_name in mappings:
        if mo_bc == mo and gid == group_index:
            return mark_name

    for gid, mo_bc, mark_name in mappings:
        if mo_bc == mo and gid == -1:
            return mark_name

    return None


def build_output_targets(group_names, mark_names, mappings):
    targets = []
    used = set()

    for group_name in group_names:
        group_index = group_names.index(group_name)
        for mark_name in mark_names:
            for gid, _, mapped_mark in mappings:
                if mapped_mark != mark_name:
                    continue
                if gid == -1 or gid == group_index:
                    key = (group_name, mark_name)
                    if key not in used:
                        targets.append(key)
                        used.add(key)
                    break

    return targets


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
    mappings, mark_names = load_mo_map(args.mo_map, args.sample, group_names)
    targets = build_output_targets(group_names, mark_names, mappings)

    r1_handles = OrderedDict()
    r2_handles = OrderedDict()
    target_barcodes = {}

    for group_name, mark_name in targets:
        stem = f"{args.sample}_{group_name}_{mark_name}"
        r1_handles[(group_name, mark_name)] = gzip.open(args.output_dir / f"{stem}_R1.fq.gz", "wt")
        r2_handles[(group_name, mark_name)] = gzip.open(args.output_dir / f"{stem}_R2.fq.gz", "wt")
        target_barcodes[(group_name, mark_name)] = set()

    try:
        for r1_rec, r2_rec in zip(fastq_iter(args.r1), fastq_iter(args.r2)):
            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, r2_comment = parse_header(r2_rec[0])
            if r1_name != r2_name:
                raise ValueError(f"Read name mismatch: {r1_name} != {r2_name}")

            if "NoMatch" in r1_comment:
                continue

            cb = find_tag_value(r1_comment, "CB")
            mo = find_tag_value(r1_comment, "MO")
            sb = find_tag_value(r1_comment, "SB")
            if not cb or not mo or not sb:
                raise ValueError(f"Missing CB or MO or SB tag in FASTQ comment for {r1_name}")

            group_name = resolve_group(args.sample, sb, sb_to_group)
            mark_name = find_mark_for_mo(mo, group_name, group_names, mappings)
            if mark_name is None:
                raise ValueError(f"MO barcode not found for sample {args.sample}: {mo}")

            key = (group_name, mark_name)
            write_fastq_record(r1_handles[key], r1_name, r1_comment, r1_rec[1], r1_rec[3])
            write_fastq_record(r2_handles[key], r2_name, r2_comment, r2_rec[1], r2_rec[3])
            target_barcodes[key].add(cb)
    finally:
        for handle in r1_handles.values():
            handle.close()
        for handle in r2_handles.values():
            handle.close()

    for group_name, mark_name in targets:
        header_path = args.output_dir / f"SAM_RG_Header_{args.sample}_{group_name}_{mark_name}.tsv"
        write_rg_header(header_path, args.sample, args.library_name, target_barcodes[(group_name, mark_name)])


def real_split(args):
    if shutil.which("codon") is None:
        raise RuntimeError("codon executable not found in PATH")

    with tempfile.TemporaryDirectory(prefix="tresflow_split_reads_dna_") as tmpdir:
        tmp_path = Path(tmpdir)
        cmd = [
            "codon",
            "run",
            "-plugin",
            "seq",
            "-release",
            str(args.script),
            args.sample,
            str(tmp_path),
            args.library_name,
            "dna",
            str(args.mo_map),
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
            raise RuntimeError(f"No DNA split outputs were produced for sample {args.sample}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["real", "mock"])
    parser.add_argument("--script", required=True, type=Path)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--mo-map", required=True, type=Path)
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
