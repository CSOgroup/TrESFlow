#!/usr/bin/env python3
"""Use the same complete primary target fragments for CPM, MACS3 and overlap QC."""
import argparse
import csv
import gzip
import json
import re
import subprocess
import sys
from bisect import bisect_right
from collections import Counter, defaultdict
from itertools import groupby
from pathlib import Path

import pyBigWig
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from tresflow_qc.core import read_duplicate_metrics


def fragments(records):
    mates = defaultdict(list)
    for read in records:
        if read.flag & 0xB0C or not read.is_proper_pair or read.reference_id != read.next_reference_id:
            continue
        if read.is_read1 == read.is_read2:
            continue
        ends = (read.reference_start, read.next_reference_start)
        key = (read.reference_id, read.get_tag("RG"), read.get_tag("CB"), read.get_tag("SB"),
               ends if read.is_read1 else ends[::-1], read.is_read1)
        mates[key].append(read)
    for key, first in mates.items():
        if key[-1]:
            yield from zip(first, mates.get((*key[:-1], False), ()))


def create_fragments(by_name, sizes):
    lengths = Counter()
    duplicates = 0
    with pysam.AlignmentFile(by_name, "rb") as bam, open("target.bedpe", "w") as bed, \
            pysam.AlignmentFile("spans.bam", "wb", template=bam) as coverage:
        for _, records in groupby(bam, key=lambda r: r.query_name):
            for r1, r2 in fragments(records):
                chrom = r1.reference_name
                left, right = sorted((r1, r2), key=lambda r: r.reference_start)
                lo, hi = left.reference_start, max(left.reference_end, right.reference_end)
                if not 0 <= lo < hi <= sizes[chrom]:
                    raise ValueError(f"{chrom}:{lo}-{hi} outside configured chromosome sizes")
                # MACS3 BEDPE is three columns describing the full fragment,
                # rather than the ten-column BEDPE format used by bedtools.
                print(chrom, lo, hi, sep="\t", file=bed)
                # deepTools gets one exact contiguous interval per eligible fragment.
                # This also covers long inserts and overlapping mates exactly once.
                span = pysam.AlignedSegment(bam.header)
                span.query_name = r1.query_name
                span.reference_id = r1.reference_id
                span.reference_start = lo
                span.mapping_quality = 60
                span.cigartuples = [(0, hi - lo)]
                coverage.write(span)
                lengths[hi - lo] += 1
                duplicates += r1.is_duplicate or r2.is_duplicate
    return lengths, duplicates


def bed_intervals(path):
    intervals = defaultdict(list)
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as handle:
        for raw in handle:
            if raw.strip() and not raw.startswith(("#", "track ", "browser ")):
                chrom, lo, hi = raw.split()[:3]
                lo, hi = int(lo), int(hi)
                if not 0 <= lo < hi:
                    raise ValueError(f"{path}: invalid BED interval {raw.strip()}")
                intervals[chrom].append((lo, hi))
    result = {}
    for chrom, spans in intervals.items():
        starts, ends = [], []
        for lo, hi in sorted(spans):
            if ends and ends[-1] >= lo:
                ends[-1] = max(ends[-1], hi)
            else:
                starts.append(lo)
                ends.append(hi)
        result[chrom] = starts, ends
    return result


def overlap_count(path):
    intervals = bed_intervals(path)
    count = 0
    with open("target.bedpe") as handle:
        for raw in handle:
            chrom, lo, hi = raw.split()
            starts, ends = intervals.get(chrom, ([], []))
            i = bisect_right(ends, int(lo))
            count += i < len(starts) and starts[i] < int(hi)
    bases = sum(sum(ends) - sum(starts) for starts, ends in intervals.values())
    return count, bases


def call_peaks(args, n, sizes, peak_bed):
    command = [args.macs3, "callpeak", "--treatment", "target.bedpe", "--format", "BEDPE",
               "--gsize", str(args.effective), "--keep-dup", "all", "--nolambda",
               "--qvalue", "0.05", "--broad", "--broad-cutoff", "0.1", "--name", "target"]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout, flush=True)
    if result.returncode:
        raise RuntimeError(f"MACS3 failed with exit status {result.returncode}")
    # The XLS header records input fragments even when no peaks are called.
    valid = re.search(r"^# total fragments in treatment: (\d+)\s*$",
                      Path("target_peaks.xls").read_text(), re.MULTILINE)
    if not valid or int(valid.group(1)) != n:
        raise ValueError(f"MACS3's effective treatment count does not reconcile to {n} usable fragments")
    count = 0
    with open("target_peaks.broadPeak") as handle, peak_bed.open("w") as out:
        for raw in handle:
            if not raw.strip():
                continue
            chrom, lo, hi = raw.split()[:3]
            start, end = int(lo), int(hi)
            # MACS3 broadPeak coordinates are already zero-based, half-open.
            if chrom not in sizes or not 0 <= start < end <= sizes[chrom]:
                raise ValueError(f"MACS3 returned an invalid peak interval: {raw.strip()}")
            print(chrom, start, end, sep="\t", file=out)
            count += 1
    return count


def analyze(args):
    sizes = dict((c, int(n)) for c, n in (line.split() for line in args.sizes.read_text().splitlines()))
    if not 0 < args.effective <= sum(sizes.values()):
        raise ValueError("effective genome size must be within the configured chromosome sizes")
    lengths, duplicates = create_fragments(args.by_name, sizes)
    n = sum(lengths.values())
    peaks, bw = Path(f"{args.split}.peaks.bed"), f"{args.split}.cpm.bw"
    if n:
        pysam.sort("-@", str(args.threads), "-o", "coverage.bam", "spans.bam")
        pysam.index("coverage.bam")
        subprocess.run([args.coverage, "-b", "coverage.bam", "-o", bw, "--normalizeUsing", "CPM",
                        "--binSize", "50", "--numberOfProcessors", str(args.threads)], check=True)
        count = call_peaks(args, n, sizes, peaks)
    else:
        with pyBigWig.open(bw, "w") as empty:
            empty.addHeader(list(sizes.items()))
        peaks.write_text("")
        count = 0
    overlap, bases = overlap_count(peaks)
    with args.retention.open() as h:
        qc = {row["metric"]: int(row["pairs"]) for row in csv.DictReader(h, delimiter="\t")}
    picard = read_duplicate_metrics(args.metrics)
    qc.update({k: picard.get(k) for k in ("read_pairs_examined", "read_pair_duplicates", "read_pair_optical_duplicates")})
    qc["picard"] = picard
    qc.update(usable_fragments=n, usable_duplicate_fragments=duplicates, peak_count=count, peak_bases=bases,
              fragments_in_peaks=overlap, frip=overlap/n if n else None, fragment_lengths=lengths,
              caller_status="completed" if n else "no usable fragments", coverage_bin_bp=50,
              effective_genome_size=args.effective, peak_caller="MACS3", macs3_format="BEDPE",
              macs3_broad=True, macs3_q_value=0.05, macs3_broad_cutoff=0.1,
              macs3_keep_duplicates=True, macs3_nolambda=True)
    if args.reference:
        qc["reference_bed"] = args.reference.name
        qc["fragments_in_reference"] = overlap_count(args.reference)[0]
    qc["depth"] = [[f * qc["bwa_primary_pairs"], f * n] for f in (0.1, 0.25, 0.5, 0.75, 1)]
    Path(f"{args.split}.dna.qc.json").write_text(json.dumps(qc, sort_keys=True))
    for name in ("spans.bam", "coverage.bam", "coverage.bam.bai", "target.bedpe",
                 "target_peaks.xls", "target_peaks.broadPeak", "target_peaks.gappedPeak"):
        Path(name).unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("split")
    for field in ("by_name", "sizes", "retention", "metrics"):
        parser.add_argument(field, type=Path)
    parser.add_argument("effective", type=int)
    parser.add_argument("threads", type=int)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--macs3", default="macs3")
    parser.add_argument("--coverage", default="bamCoverage")
    analyze(parser.parse_args())
