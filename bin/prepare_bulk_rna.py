#!/usr/bin/env python3
"""Validate and prepare one routed RNA pair stream directly as STAR input."""
import base64
import json
import sys
from collections import Counter
from itertools import zip_longest
from pathlib import Path

from export_rna_fastqs import pool_identifier, retained_sequences
from preprocessed_fastqs import CLIPPED_RNA_SUFFIX, observe_routed, read_tags, validate_pair
from tresflow_fastq_utils import fastq_iter, normalize_qname, parse_header


def prepare(spec, r1, r2, output):
    counts = {field: Counter() for field in ("L1", "L2", "L3", "SB")}
    total = kept = 0
    width = spec["ligation_length"] // 3
    prefix = max(spec["sb_start"] + spec["sb_length"], spec["umi_start"] + spec["umi_length"])
    with output.open("w") as sam:
        sam.write("@HD\tVN:1.6\tSO:unsorted\n")
        for total, (first, second) in enumerate(zip_longest(fastq_iter(r1), fastq_iter(r2)), 1):
            try:
                if first is None or second is None:
                    raise ValueError("unequal mate counts")
                validate_pair(first, second, spec, spec.get("sb_indices", {}), width * 3)
                name, comment = parse_header(first[0])
                name = normalize_qname(name)
                tags = read_tags(comment)
                suffix = CLIPPED_RNA_SUFFIX.search(name)
                if suffix:
                    pool = suffix[1]
                    name = name[:suffix.start()]
                else:
                    if "SB" not in tags:
                        # The physical oligo can be remapped to any logical index.
                        # Only an exact-match tag layout can recover it from R2.
                        start, length = spec["sb_start"], spec["sb_length"]
                        bases = second[1][start:start + length].upper()
                        if spec["sb_hd"] or len(bases) != length or set(bases) - set("ACGT"):
                            raise ValueError("missing physical SB tag and ambiguous R2 sample-barcode prefix")
                        tags["SB"] = bases.translate(str.maketrans("ACGT", "TGCA"))[::-1] if spec["sb_rc"] else bases
                    pool = pool_identifier(tags["CB"], tags["SB"])
                observe_routed(counts, tags, width)
                pair = retained_sequences(first, second, 0 if suffix else prefix, 20)
                if pair is None:
                    continue
                kept += 1
                # STAR writes CB=XP and UB after correction. XI preserves the full TrES identity.
                attrs = [f"{tag}:Z:{tags[tag]}" for tag in ("XI", "SB", "UM") if tag in tags]
                attrs.extend([f"XP:Z:{pool}", f"ZI:i:{kept}"])
                for flag, (seq, qual) in zip((77, 141), pair):
                    print(name, flag, "*", 0, 0, "*", "*", 0, 0, seq, qual, *attrs, sep="\t", file=sam)
            except ValueError as exc:
                raise ValueError(f'{spec.get("r1", r1)}, {spec.get("r2", r2)}, pair {total}: {exc}') from exc
    if spec.get("pairs") is not None and total != spec["pairs"]:
        raise ValueError(f"{r1}: expected {spec['pairs']} pairs, found {total}")
    return dict(routed_pairs=total, retained_pairs=kept, barcodes=counts)


if __name__ == "__main__":
    spec = json.loads(base64.b64decode(sys.argv[1]))
    qc = prepare(spec, *map(Path, sys.argv[2:5]))
    Path(sys.argv[5]).write_text(json.dumps(qc))
