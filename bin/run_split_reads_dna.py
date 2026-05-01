#!/usr/bin/env python3

import argparse
import subprocess
import sys
import tempfile
import time
from collections import OrderedDict
from pathlib import Path

from tresflow_fastq_utils import (
    canonicalize_fastq_comment,
    compress_fastq_with_pigz,
    compress_final_fastqs,
    fastq_iter,
    find_tag_value,
    load_sb_group_map,
    log_event,
    move_split_output,
    normalize_split_fastq_name,
    parse_header,
    resolve_codon_bin,
    resolve_group,
    resolve_temp_root,
    write_fastq_record,
    write_rg_header,
)


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


def mock_split(args):
    sb_to_group, group_names = load_sb_group_map(args.sb_group_map, args.sample)
    mappings, mark_names = load_mo_map(args.mo_map, args.sample, group_names)
    targets = build_output_targets(group_names, mark_names, mappings)

    r1_handles = OrderedDict()
    r2_handles = OrderedDict()
    target_barcodes = {}

    for group_name, mark_name in targets:
        stem = f"{args.sample}_{group_name}_{mark_name}"
        r1_handles[(group_name, mark_name)] = open(
            args.output_dir / f"{stem}_R1.fastq",
            "wt",
            encoding="utf-8",
        )
        r2_handles[(group_name, mark_name)] = open(
            args.output_dir / f"{stem}_R2.fastq",
            "wt",
            encoding="utf-8",
        )
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
            r1_comment = canonicalize_fastq_comment(args.sample, group_name, r1_comment)
            r2_comment = canonicalize_fastq_comment(args.sample, group_name, r2_comment)
            canonical_cb = find_tag_value(r1_comment, "CB")
            write_fastq_record(r1_handles[key], r1_name, r1_comment, r1_rec[1], r1_rec[3])
            write_fastq_record(r2_handles[key], r2_name, r2_comment, r2_rec[1], r2_rec[3])
            target_barcodes[key].add(canonical_cb)
    finally:
        for handle in r1_handles.values():
            handle.close()
        for handle in r2_handles.values():
            handle.close()

    for group_name, mark_name in targets:
        header_path = args.output_dir / f"SAM_RG_Header_{args.sample}_{group_name}_{mark_name}.tsv"
        write_rg_header(header_path, args.sample, args.library_name, target_barcodes[(group_name, mark_name)])

    compress_final_fastqs(args.output_dir, args.pigz_threads)


def real_split(args):
    codon_bin = resolve_codon_bin()

    with tempfile.TemporaryDirectory(prefix="tresflow_split_reads_dna_", dir=resolve_temp_root()) as tmpdir:
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
            "dna",
            str(args.mo_map),
            str(args.r1),
            str(args.r2),
            str(args.sb_group_map),
        ]
        codon_start = time.monotonic()
        log_event("Starting Codon Split_ReadsV2.codon DNA", args.r1, args.r2)
        subprocess.run(cmd, check=True)
        log_event("Finished Codon Split_ReadsV2.codon DNA", args.r1, args.r2, elapsed=time.monotonic() - codon_start)

        moved = 0
        fastq_moved = 0
        for pattern in (
            f"{args.sample}_*_R1.fastq.gz",
            f"{args.sample}_*_R2.fastq.gz",
            f"{args.sample}_*_R1.fq.gz",
            f"{args.sample}_*_R2.fq.gz",
            f"{args.sample}_*_R1.fastq",
            f"{args.sample}_*_R2.fastq",
            f"{args.sample}_*_R1.fq",
            f"{args.sample}_*_R2.fq",
            f"SAM_RG_Header_{args.sample}_*.tsv",
        ):
            for source in sorted(tmp_path.glob(pattern)):
                destination = move_split_output(source, args.output_dir)
                moved += 1
                if destination.name.endswith((".fastq.gz", ".fastq")):
                    fastq_moved += 1

        if moved == 0 or fastq_moved == 0:
            raise RuntimeError(f"No DNA split outputs were produced for sample {args.sample}")

        compress_final_fastqs(args.output_dir, args.pigz_threads)


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
    parser.add_argument("--pigz-threads", required=True, type=int)
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
