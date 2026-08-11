#!/usr/bin/env python3
"""Count nested RNA BAM filter populations from one SAM stream.

One primary R1 alignment is used as the stable representative of each input
read pair. This matches the pair unit reported as `read1` by samtools flagstat
while avoiding double-counting STAR secondary/supplementary alignments.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_values(path: Path) -> set[str]:
    values = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            value = line.strip().split("\t", 1)[0]
            if value and not value.startswith("#"):
                values.add(value)
    return values


def rg_value(fields: list[str]) -> str | None:
    for field in fields[11:]:
        if field.startswith("RG:Z:"):
            return field[5:]
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-id", required=True)
    parser.add_argument("--canonical-contigs", required=True, type=Path)
    parser.add_argument("--called-barcodes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    canonical = read_values(args.canonical_contigs)
    called = read_values(args.called_barcodes)
    counts = {
        "star_mapped_primary_pairs": 0,
        "paired_filter_pairs": 0,
        "canonical_pairs": 0,
        "called_cell_pairs": 0,
    }

    for raw_line in sys.stdin:
        if not raw_line or raw_line.startswith("@"):
            continue
        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 11:
            raise ValueError("Malformed SAM record in RNA retention stream")
        flag = int(fields[1])

        # Fixed read-pair representative: primary, non-supplementary R1.
        if not (flag & 0x40) or (flag & 0x100) or (flag & 0x800):
            continue
        if (flag & 0x4) or fields[2] == "*":
            continue
        counts["star_mapped_primary_pairs"] += 1

        # The current pipeline passes `0x1,0x2` to samtools. Samtools parses
        # that numeric comma expression as 0x1, so the effective existing
        # criterion is paired, not proper-pair. This script reports, but does
        # not alter, that behavior.
        if not (flag & 0x1):
            continue
        counts["paired_filter_pairs"] += 1

        if fields[2] not in canonical:
            continue
        counts["canonical_pairs"] += 1

        if rg_value(fields) not in called:
            continue
        counts["called_cell_pairs"] += 1

    ordered = [
        "star_mapped_primary_pairs",
        "paired_filter_pairs",
        "canonical_pairs",
        "called_cell_pairs",
    ]
    for previous, current in zip(ordered, ordered[1:]):
        if counts[current] > counts[previous]:
            raise ValueError(
                f"Non-nested RNA retention counts: {current}={counts[current]} "
                f"> {previous}={counts[previous]}"
            )

    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("split_id\tmetric\tpairs\tunit\n")
        for metric in ordered:
            handle.write(
                f"{args.split_id}\t{metric}\t{counts[metric]}\t"
                "primary_read1_pair_representatives\n"
            )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        sys.stderr.write(f"ERROR: {error}\n")
        raise SystemExit(1)
