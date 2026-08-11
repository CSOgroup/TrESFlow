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
    fastq_iter,
    find_tag_value,
    load_sb_group_map,
    log_event,
    move_split_output,
    parse_header,
    resolve_codon_bin,
    resolve_group,
    resolve_temp_root,
    write_fastq_record,
    write_rg_header,
)


def mock_split(args):
    sb_to_group, group_names = load_sb_group_map(args.sb_group_map, args.sample)

    r1_handles = OrderedDict()
    r2_handles = OrderedDict()
    group_barcodes = {}

    for group_name in group_names:
        r1_handles[group_name] = open(
            args.output_dir / f"{args.sample}_{group_name}_R1.fastq",
            "wt",
            encoding="utf-8",
        )
        r2_handles[group_name] = open(
            args.output_dir / f"{args.sample}_{group_name}_R2.fastq",
            "wt",
            encoding="utf-8",
        )
        group_barcodes[group_name] = set()

    processed = 0
    accepted = 0
    group_counts = OrderedDict((group_name, 0) for group_name in group_names)

    try:
        for r1_rec, r2_rec in zip(fastq_iter(args.r1), fastq_iter(args.r2)):
            processed += 1
            r1_name, r1_comment = parse_header(r1_rec[0])
            r2_name, r2_comment = parse_header(r2_rec[0])
            if r1_name != r2_name:
                raise ValueError(f"Read name mismatch: {r1_name} != {r2_name}")

            if "NoMatch" in r1_comment:
                continue
            accepted += 1

            cb = find_tag_value(r1_comment, "CB")
            sb = find_tag_value(r1_comment, "SB")
            if not cb or not sb:
                raise ValueError(f"Missing CB or SB tag in FASTQ comment for {r1_name}")

            group_name = resolve_group(args.sample, sb, sb_to_group)
            group_counts[group_name] += 1
            r1_comment = canonicalize_fastq_comment(args.sample, group_name, r1_comment)
            r2_comment = canonicalize_fastq_comment(args.sample, group_name, r2_comment)
            canonical_cb = find_tag_value(r1_comment, "CB")
            write_fastq_record(r1_handles[group_name], r1_name, r1_comment, r1_rec[1], r1_rec[3])
            write_fastq_record(r2_handles[group_name], r2_name, r2_comment, r2_rec[1], r2_rec[3])
            group_barcodes[group_name].add(canonical_cb)
    finally:
        for handle in r1_handles.values():
            handle.close()
        for handle in r2_handles.values():
            handle.close()

    for group_name in group_names:
        header_path = args.output_dir / f"SAM_RG_Header_{args.sample}_{group_name}.tsv"
        write_rg_header(header_path, args.sample, args.library_name, group_barcodes[group_name])

    metrics_path = args.output_dir / f"{args.sample}.rna_read_retention.tsv"
    with metrics_path.open("wt", encoding="utf-8") as handle:
        handle.write("sample_id\tmodality\tgroup\tbranch\tmetric\tpairs\tunit\n")
        handle.write(f"{args.sample}\trna\t__all__\t__all__\ttrimmed_input_pairs\t{processed}\tread_pairs\n")
        handle.write(f"{args.sample}\trna\t__all__\t__all__\tjoint_barcode_accepted_pairs\t{accepted}\tread_pairs\n")
        for group_name, count in group_counts.items():
            handle.write(f"{args.sample}\trna\t{group_name}\t__all__\trouted_group_pairs\t{count}\tread_pairs\n")
            handle.write(f"{args.sample}\trna\t{group_name}\tRNA\trouted_branch_pairs\t{count}\tread_pairs\n")


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
        codon_start = time.monotonic()
        log_event("Starting Codon Split_ReadsV2.codon RNA", args.r1, args.r2)
        subprocess.run(cmd, check=True)
        log_event("Finished Codon Split_ReadsV2.codon RNA", args.r1, args.r2, elapsed=time.monotonic() - codon_start)

        retention_name = f"{args.sample}.rna_read_retention.tsv"
        retention_source = tmp_path / retention_name
        if not retention_source.is_file():
            raise FileNotFoundError(
                f"Codon completed without producing required RNA retention metrics: {retention_source}"
            )

        moved = 0
        fastq_moved = 0
        for pattern in (
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
            raise RuntimeError(f"No RNA split outputs were produced for sample {args.sample}")

        retention_destination = move_split_output(retention_source, args.output_dir)
        if retention_destination != args.output_dir / retention_name:
            raise RuntimeError(
                f"RNA retention metrics were moved to an unexpected path: {retention_destination}"
            )


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
