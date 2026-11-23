#!/usr/bin/env python3
"""Unstranded coverage of each representative template's union of aligned blocks."""
import json
import sys
from collections import defaultdict
from itertools import groupby
from pathlib import Path

import pyBigWig
import pysam


def blocks():
    with pysam.AlignmentFile("-", "rb") as bam:
        for _, group in groupby(bam, key=lambda r: r.query_name):
            spans = defaultdict(list)
            for read in group:
                if not read.is_unmapped:
                    spans[read.get_tag("ZI"), read.reference_id].extend(read.get_blocks())
            for (_, chrom), intervals in spans.items():
                merged = []
                for start, end in sorted(intervals):
                    if merged and start <= merged[-1][1]:
                        merged[-1] = merged[-1][0], max(end, merged[-1][1])
                    else:
                        merged.append((start, end))
                for lo, hi in merged:
                    print(chrom, lo, 1, sep="\t")
                    print(chrom, hi, -1, sep="\t")


def coverage(bam_path, qc, output):
    n = json.loads(Path(qc).read_text())["molecules"]
    with pysam.AlignmentFile(bam_path, "rb") as bam, pyBigWig.open(output, "w") as track:
        track.addHeader(list(zip(bam.references, bam.lengths)))
        events = (tuple(map(int, line.split())) for line in sys.stdin)
        for chrom, chromosome in groupby(events, key=lambda row: row[0]):
            depth = previous = 0
            starts, ends, values = [], [], []

            def flush():
                if starts:
                    track.addEntries([bam.references[chrom]] * len(starts), starts, ends=ends, values=values)
                    starts.clear()
                    ends.clear()
                    values.clear()

            for pos, at_position in groupby(chromosome, key=lambda row: row[1]):
                if depth:
                    starts.append(previous)
                    ends.append(pos)
                    values.append(depth * 1e6 / n)
                    if len(starts) == 8192:
                        flush()
                depth += sum(row[2] for row in at_position)
                previous = pos
            flush()
            if depth:
                raise ValueError(f"Unbalanced coverage events on {bam.references[chrom]}")


if __name__ == "__main__":
    if sys.argv[1] == "blocks":
        blocks()
    else:
        coverage(*sys.argv[2:])
