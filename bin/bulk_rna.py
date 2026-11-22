#!/usr/bin/env python3
"""Reconcile primary paired representatives with STAR's unique GeneFull matrix."""
import csv
import json
import sqlite3
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from itertools import groupby
from pathlib import Path

import pysam


class Regions:
    """Merged exon and gene intervals in STAR's concatenated reference coordinates."""
    def __init__(self, index):
        names = (index / "chrName.txt").read_text().splitlines()
        self.offset = dict(zip(names, map(int, (index / "chrStart.txt").read_text().split())))
        intervals = defaultdict(list)
        genes = {}
        with (index / "exonGeTrInfo.tab").open() as handle:
            next(handle)
            for line in handle:
                start, end, strand, gene, _ = map(int, line.split())
                intervals["exon", strand].append((start, end + 1))
                lo, hi = genes.get((gene, strand), (start, end + 1))
                genes[gene, strand] = min(lo, start), max(hi, end + 1)
        for (_, strand), span in genes.items():
            intervals["gene", strand].append(span)
        self.bounds = {}
        for category, spans in intervals.items():
            starts, ends = [], []
            for lo, hi in sorted(set(spans)):
                if ends and lo <= ends[-1]:
                    ends[-1] = max(ends[-1], hi)
                else:
                    starts.append(lo)
                    ends.append(hi)
            self.bounds[category] = starts, ends

    def classify(self, reads):
        blocks = [(self.offset[r.reference_name] + lo, self.offset[r.reference_name] + hi)
                  for r in reads for lo, hi in r.get_blocks()]
        for kind, label in (("exon", "exonic"), ("gene", "intronic")):
            strands = set()
            for strand in (0, 1, 2):
                starts, ends = self.bounds.get((kind, strand), ([], []))
                if any(bisect_left(starts, hi) > bisect_right(ends, lo) for lo, hi in blocks):
                    strands.add(strand)
            if strands:
                orientation = "ambiguous"
                if len(strands) == 1 and 0 not in strands:
                    read = reads[0]
                    forward = read.is_reverse != read.is_read1
                    orientation = "sense" if (1 if forward else 2) in strands else "antisense"
                return label, orientation
        return "intergenic", None


def matrix(raw):
    genes = [line.split("\t", 1)[0] for line in (raw / "features.tsv").read_text().splitlines()]
    pools = (raw / "barcodes.tsv").read_text().splitlines()
    expected = Counter()
    with (raw / "matrix.mtx").open() as handle:
        data = (line for line in handle if not line.startswith("%"))
        n, m, entries = map(int, next(data).split())
        if (n, m) != (len(genes), len(pools)):
            raise ValueError("STAR matrix dimensions disagree with features and barcodes")
        for line in data:
            gene, pool, count = map(int, line.split())
            expected[genes[gene - 1], pools[pool - 1]] = count
        if len(expected) != entries:
            raise ValueError("incomplete STAR matrix")
    return genes, expected


def summarize(split, index, prep, source, solo):
    genes, expected = matrix(solo / "raw")
    annotation = Regions(index)
    stats = Counter()
    qc = json.loads(prep.read_text())
    database = sqlite3.connect("molecules.sqlite")
    database.execute("CREATE TABLE molecules (gene, pool, umi, reads INTEGER, PRIMARY KEY(gene,pool,umi)) WITHOUT ROWID")
    actual = Counter()
    # Name-sort order defines the representative. ZI distinguishes repeated
    # sequencing names in different technical chunks of the same library.
    with pysam.AlignmentFile(source, "rb") as bam, pysam.AlignmentFile("representatives.bam", "wb", template=bam) as out:
        for _, records in groupby(bam, key=lambda r: r.query_name):
            by_input = defaultdict(list)
            for r in records:
                if not (r.is_unmapped or r.is_secondary or r.is_supplementary):
                    by_input[r.get_tag("ZI")].append(r)
            for _, reads in sorted(by_input.items()):
                one = reads[0]
                if one.get_tag("NH") != 1:
                    continue
                region, strand = annotation.classify(reads)
                stats[region] += 1
                if strand:
                    stats[strand] += 1
                if not all(one.has_tag(tag) and one.get_tag(tag) not in ("", "-") for tag in ("GX", "UB", "CB")):
                    continue
                key = tuple(one.get_tag(tag) for tag in ("GX", "CB", "UB"))
                if key[:2] not in expected:
                    continue
                if len(reads) != 2 or {r.is_read1 for r in reads} != {True, False}:
                    raise ValueError(f"{one.query_name}: counted GeneFull molecule lacks a complete primary pair")
                if any(tuple(r.get_tag(t) for t in ("GX", "XP", "UB")) != key for r in reads):
                    raise ValueError(f"{one.query_name}: inconsistent STAR gene/pool/corrected UMI")
                inserted = database.execute("INSERT OR IGNORE INTO molecules VALUES (?,?,?,1)", key).rowcount
                if inserted:
                    actual[key[:2]] += 1
                    for r in reads:
                        r.set_tag("ZI", None)
                        out.write(r)
                else:
                    database.execute("UPDATE molecules SET reads=reads+1 WHERE gene=? AND pool=? AND umi=?", key)
    if actual != expected:
        mismatch = next(key for key in expected.keys() | actual.keys() if expected[key] != actual[key])
        raise ValueError(f"GeneFull {mismatch}: matrix={expected[mismatch]}, BAM molecules={actual[mismatch]}")

    counts = Counter()
    for (gene, _), count in actual.items():
        counts[gene] += count
    with open(f"{split}.counts.tsv", "w") as handle:
        print("gene_id", split, sep="\t", file=handle)
        for gene in genes:
            print(gene, counts[gene], sep="\t", file=handle)

    def yield_at(group, fraction):
        return sum(n * (1 - (1 - fraction) ** k) for k, n in group)

    molecule_support = list(database.execute("SELECT reads,COUNT(*) FROM molecules GROUP BY reads"))
    gene_support = Counter(n for (n,) in database.execute("SELECT SUM(reads) FROM molecules GROUP BY gene"))
    database.close()
    qc["depth"] = [[f * qc["retained_pairs"], yield_at(molecule_support, f), yield_at(gene_support.items(), f)]
                   for f in (0.1, 0.25, 0.5, 0.75, 1)]
    qc.update(stats, molecules=sum(counts.values()), genes_detected=len(counts), pools_detected=len({p for _, p in actual}))
    with (solo / "CellReads.stats").open() as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    for field in ("genomeU", "genomeM", "featureU", "countedU"):
        qc[field] = sum(int(row[field]) for row in rows)
    Path(f"{split}.rna.qc.json").write_text(json.dumps(qc, sort_keys=True))


if __name__ == "__main__":
    summarize(sys.argv[1], *map(Path, sys.argv[2:]))
