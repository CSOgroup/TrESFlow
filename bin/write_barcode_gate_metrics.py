#!/usr/bin/env python3
"""Write exact cumulative barcode-gate and composition counters.

Barcode matching happens before paired trimming in TrESFlow.  The tag-record
stream records those existing decisions for every raw pair.  This metrics-only
helper intersects that stream with the paired FASTQs which actually reach
Split_ReadsV2.  It never repeats barcode matching and never writes FASTQ data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter, OrderedDict
from itertools import zip_longest
from pathlib import Path

from run_split_reads_dna import find_mark_for_mo, load_mo_map
from tresflow_fastq_utils import (
    fastq_iter,
    load_sb_group_map,
    parse_header,
    resolve_group,
)


NO_MATCH = "NoMatch"


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("rt", encoding="utf-8")


def tag_value(fields: list[str], tag: str) -> str:
    prefix = f"{tag}:Z:"
    for field in fields:
        if field.startswith(prefix):
            return field[len(prefix) :]
    raise ValueError(f"Tag record is missing required {tag}:Z field")


def tag_record_iter(path: Path):
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                raise ValueError(f"Malformed tag record at {path}:{line_number}")
            yield fields[0], fields[1:]


def retained_qnames(r1: Path, r2: Path):
    for pair_number, pair in enumerate(zip_longest(fastq_iter(r1), fastq_iter(r2)), start=1):
        r1_record, r2_record = pair
        if r1_record is None or r2_record is None:
            raise ValueError(f"Paired FASTQs have different record counts at pair {pair_number}")
        r1_name, _ = parse_header(r1_record[0])
        r2_name, _ = parse_header(r2_record[0])
        if r1_name != r2_name:
            raise ValueError(f"Paired FASTQ name mismatch at pair {pair_number}: {r1_name} != {r2_name}")
        yield r1_name


def matched_tag_fields(qnames, records):
    record_iterator = iter(records)
    current = next(record_iterator, None)
    for qname in qnames:
        while current is not None and current[0] != qname:
            current = next(record_iterator, None)
        if current is None:
            raise ValueError(f"Retained FASTQ QNAME was not found in ordered tag records: {qname}")
        yield qname, current[1]
        current = next(record_iterator, None)


def percentage(count: int, denominator: int) -> float:
    return count / denominator * 100.0 if denominator else 0.0


def write_gate_rows(path: Path, sample: str, modality: str, counts: OrderedDict[str, int]) -> None:
    denominator = counts["split_input_pairs"]
    definitions = {
        "split_input_pairs": "paired reads surviving all preceding FASTQ processing and entering splitting/routing",
        "ligation_barcode_accepted_pairs": "split-input pairs with L1, L2 and L3 accepted on the same pair",
        "sample_barcode_accepted_pairs": "ligation-accepted pairs with an accepted configured sample barcode",
        "modality_barcode_accepted_pairs": "sample-barcode-accepted DNA pairs with an accepted configured MO barcode",
    }
    with path.open("wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "modality",
                "metric",
                "pairs",
                "unit",
                "denominator_metric",
                "denominator_pairs",
                "definition",
                "source",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for metric, count in counts.items():
            writer.writerow(
                {
                    "sample_id": sample,
                    "modality": modality,
                    "metric": metric,
                    "pairs": count,
                    "unit": "read_pairs",
                    "denominator_metric": "split_input_pairs",
                    "denominator_pairs": denominator,
                    "definition": definitions[metric],
                    "source": "existing_tag_decisions_intersected_with_split_input_qnames",
                }
            )


def write_composition_rows(
    path: Path,
    sample: str,
    modality: str,
    sample_counts: OrderedDict[str, int],
    sample_denominator: int,
    mark_counts: dict[str, OrderedDict[str, int]],
    mark_denominators: dict[str, int],
) -> None:
    fields = [
        "sample_id",
        "modality",
        "group",
        "barcode_type",
        "category",
        "count",
        "percentage",
        "denominator_count",
        "denominator_definition",
        "source",
    ]
    with path.open("wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for category, count in sample_counts.items():
            writer.writerow(
                {
                    "sample_id": sample,
                    "modality": modality,
                    "group": "__all__",
                    "barcode_type": "sample_barcode",
                    "category": category,
                    "count": count,
                    "percentage": f"{percentage(count, sample_denominator):.6f}",
                    "denominator_count": sample_denominator,
                    "denominator_definition": "split-input pairs with L1, L2 and L3 accepted on the same pair",
                    "source": "existing_tag_decisions_intersected_with_split_input_qnames",
                }
            )
        for group, categories in mark_counts.items():
            denominator = mark_denominators[group]
            for category, count in categories.items():
                writer.writerow(
                    {
                        "sample_id": sample,
                        "modality": modality,
                        "group": group,
                        "barcode_type": "dna_mark",
                        "category": category,
                        "count": count,
                        "percentage": f"{percentage(count, denominator):.6f}",
                        "denominator_count": denominator,
                        "denominator_definition": (
                            f"DNA pairs assigned to sample-barcode group {group} after same-pair L1/L2/L3 acceptance"
                        ),
                        "source": "existing_tag_decisions_intersected_with_split_input_qnames",
                    }
                )


def calculate(args):
    sb_to_group, group_names = load_sb_group_map(args.sb_group_map, args.sample)
    sample_counts = OrderedDict((group, 0) for group in group_names)
    sample_counts[NO_MATCH] = 0

    mappings = None
    mark_counts: dict[str, OrderedDict[str, int]] = {}
    mark_denominators = OrderedDict((group, 0) for group in group_names)
    if args.modality == "dna":
        if args.mo_map is None:
            raise ValueError("DNA metrics require --mo-map")
        mappings, mark_names = load_mo_map(args.mo_map, args.sample, group_names)
        for group in group_names:
            configured = []
            group_index = group_names.index(group)
            for gid, _barcode, mark in mappings:
                if gid in {-1, group_index} and mark not in configured:
                    configured.append(mark)
            mark_counts[group] = OrderedDict((mark, 0) for mark in configured)
            mark_counts[group][NO_MATCH] = 0

    total = ligation = sample_accepted = modality_accepted = 0
    qnames = retained_qnames(args.r1, args.r2)
    records = tag_record_iter(args.tag_records)
    for qname, fields in matched_tag_fields(qnames, records):
        total += 1
        l_values = [tag_value(fields, f"L{layer}") for layer in (1, 2, 3)]
        if any(value == NO_MATCH for value in l_values):
            continue
        ligation += 1

        sb = tag_value(fields, "SB")
        group = None
        if sb != NO_MATCH:
            try:
                group = resolve_group(args.sample, sb, sb_to_group)
            except (KeyError, ValueError):
                group = None
        if group is None:
            sample_counts[NO_MATCH] += 1
            continue
        sample_counts[group] += 1
        sample_accepted += 1

        if args.modality == "dna":
            mark_denominators[group] += 1
            mo = tag_value(fields, "MO")
            mark = None if mo == NO_MATCH else find_mark_for_mo(mo, group, group_names, mappings)
            if mark is None:
                mark_counts[group][NO_MATCH] += 1
            else:
                mark_counts[group][mark] += 1
                modality_accepted += 1

    if sum(sample_counts.values()) != ligation:
        raise ValueError("Sample-barcode composition does not reconcile to the L1/L2/L3 denominator")
    if args.modality == "dna":
        for group in group_names:
            if sum(mark_counts[group].values()) != mark_denominators[group]:
                raise ValueError(f"DNA mark composition does not reconcile for group {group}")
        if sum(mark_denominators.values()) != sample_accepted:
            raise ValueError("DNA group denominators do not reconcile to the sample-barcode gate")

    counts = OrderedDict(
        [
            ("split_input_pairs", total),
            ("ligation_barcode_accepted_pairs", ligation),
            ("sample_barcode_accepted_pairs", sample_accepted),
        ]
    )
    if args.modality == "dna":
        counts["modality_barcode_accepted_pairs"] = modality_accepted
    values = list(counts.values())
    if any(later > earlier for earlier, later in zip(values, values[1:])):
        raise ValueError(f"Cumulative barcode gates are not nested: {counts}")
    return counts, sample_counts, mark_counts, mark_denominators


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--modality", required=True, choices=["rna", "dna"])
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--tag-records", required=True, type=Path)
    parser.add_argument("--sb-group-map", required=True, type=Path)
    parser.add_argument("--mo-map", type=Path)
    parser.add_argument("--output-gates", required=True, type=Path)
    parser.add_argument("--output-composition", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    counts, sample_counts, mark_counts, mark_denominators = calculate(args)
    write_gate_rows(args.output_gates, args.sample, args.modality, counts)
    write_composition_rows(
        args.output_composition,
        args.sample,
        args.modality,
        sample_counts,
        counts["ligation_barcode_accepted_pairs"],
        mark_counts,
        mark_denominators,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        raise SystemExit(1)
