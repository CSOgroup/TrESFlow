#!/usr/bin/env python3

import argparse
import base64
import csv
import json
import re
import shutil
import sys
from itertools import zip_longest
from pathlib import Path

from export_rna_fastqs import pool_identifier
from run_fq_to_sam import extract_rna_tags, normalize_qname, parse_header
from tresflow_fastq_utils import fastq_iter, parse_aviti_qname, write_rg_header


CLIPPED_RNA_SUFFIX = re.compile(r"_(P[0-9a-f]{64})_([^_]+)$")
LIGATION_BARCODE = re.compile(r"[ACGT]+")
UMI_BASES = re.compile(r"[ACGTN]+")
UMI_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
REQUIRED_TAGS = {"CB", "XI", "SB", "RG", "UM", "MO"}
EXACT_RNA_ADVICE = "Supply the exact tagged RNA split FASTQs from restart/rna/."


def read_tags(comment):
    tags = {}
    for token in comment.split():
        key, _, tag = token.partition(":")
        if key in REQUIRED_TAGS:
            kind, _, value = tag.partition(":")
            if kind != "Z" or not value or value == "NoMatch" or key in tags:
                raise ValueError(f"invalid or duplicate {key} tag")
            tags[key] = value
    return tags


def validate_pair(first, second, spec, barcodes, ligation_length):
    name1, comment1 = parse_header(first[0])
    name2, comment2 = parse_header(second[0])
    name = normalize_qname(name1, name2)
    if not name:
        raise ValueError(f"mate read names differ ({name1!r}, {name2!r})")
    dna = spec["modality"] == "dna"
    fastq_only = spec.get("fastq_only", False)
    clipped = None if dna else CLIPPED_RNA_SUFFIX.search(name)
    if clipped and not fastq_only:
        raise ValueError("RNA bulk export detected. Supply the original, unclipped post-split RNA FASTQs.")

    required = {"CB", "XI", "UM"}
    if dna:
        required = {"CB", "XI", "SB", "RG", "MO"}
    elif not fastq_only:
        required.update({"SB", "RG"})
    mate_tags = []
    for mate, comment in (("R1", comment1), ("R2", comment2)):
        try:
            tags = read_tags(comment)
            missing = required - tags.keys()
            if missing:
                raise ValueError(f"missing required FASTQ tags: {', '.join(sorted(missing))}")
        except ValueError as error:
            advice = f". {EXACT_RNA_ADVICE}" if not dna else ""
            raise ValueError(f"{mate}: {error}{advice}") from error
        mate_tags.append(tags)
    tags = mate_tags[0]
    if tags != mate_tags[1]:
        raise ValueError("mate FASTQ tags differ")

    namespace = f'{spec["sample"]}_{spec["group"]}_'
    cell = tags["CB"]
    sb_index, _, ligation = cell.removeprefix(namespace).partition("_")
    if (not cell.startswith(namespace) or cell != tags["XI"] or not re.fullmatch(r"[0-9]*[1-9][0-9]*", sb_index)
            or len(ligation) != ligation_length or not LIGATION_BARCODE.fullmatch(ligation)):
        advice = f". {EXACT_RNA_ADVICE}" if not dna else ""
        raise ValueError(f"CB/XI must carry the same full cell ID for samplesheet sample {spec['sample']!r} "
                         f"and group {spec['group']!r}; found CB={cell!r}, XI={tags['XI']!r}{advice}")
    if not fastq_only and barcodes.get(tags["SB"]) != sb_index:
        raise ValueError("CB/XI and corrected SB do not match the samplesheet group and logical barcode")

    if dna:
        # The Codon splitter preserves the lane token, including any leading zeros.
        parse_aviti_qname(name)
        instrument, run, flowcell, lane, *_ = name.split(":")
        expected_rg = f"{instrument}:{run}:{flowcell}:L{lane}"
        if tags["RG"] != expected_rg:
            raise ValueError(f"RG {tags['RG']!r} does not match AVITI QNAME read group {expected_rg!r}")
        mo_barcode = spec.get("mo_barcode")
        if (not fastq_only or mo_barcode is not None) and tags["MO"] != mo_barcode:
            raise ValueError(f"MO {tags['MO']!r} does not match the DNA mark barcode {mo_barcode!r}")
        return expected_rg

    if tags.get("RG", cell) != cell:
        raise ValueError("RNA RG must equal CB/XI")
    umi = tags["UM"]
    if len(umi) != spec["umi_length"] or not UMI_BASES.fullmatch(umi.upper()):
        raise ValueError(f"RNA UM must contain A/C/G/T/N and match the configured UMI length. {EXACT_RNA_ADVICE}")
    if extract_rna_tags(comment1) != extract_rna_tags(comment2):
        raise ValueError("mate RNA SAM tags differ")
    if clipped:
        pool, name_umi = clipped.groups()
        if name_umi != umi or ("SB" in tags and pool != pool_identifier(cell, tags["SB"])):
            raise ValueError("Clipped RNA QNAME pool/UMI does not match the retained FASTQ tags")
        return None

    umi_start = spec["umi_start"]
    umi_end = umi_start + spec["umi_length"]
    if fastq_only and len(second[1]) < umi_end:
        raise ValueError(f"RNA R2 is too short for the unclipped UMI layout. {EXACT_RNA_ADVICE}")
    if len(second[1]) >= umi_end:
        bases = second[1][umi_start:umi_end].upper()
        if UMI_BASES.fullmatch(bases) and bases.translate(UMI_COMPLEMENT)[::-1] != umi.upper():
            raise ValueError(f"RNA R2 prefix does not encode the UM tag. {EXACT_RNA_ADVICE}")
    return None


def validate_split(r1, r2, spec, input_rg=None, output_rg=None):
    count = 0
    groups = set()
    ligation_length = spec["ligation_length"]
    barcodes = spec.get("sb_indices", {})
    try:
        for first, second in zip_longest(fastq_iter(r1), fastq_iter(r2)):
            if first is None or second is None:
                raise ValueError("mate FASTQs have unequal record counts")
            rg = validate_pair(first, second, spec, barcodes, ligation_length)
            if rg is not None:
                groups.add(rg)
            count += 1
    except (ValueError, OSError, EOFError) as error:
        raise ValueError(f"at pair {count + 1}: {error}") from error

    expected_count = spec.get("pairs")
    if expected_count is not None and expected_count != count:
        raise ValueError(f"manifest expects {expected_count} pairs, found {count}")
    if spec["modality"] == "dna":
        write_rg_header(output_rg, spec["sample"], spec["library"], groups)
        if input_rg is not None:
            expected = {parse_rg_line(line) for line in output_rg.read_text().splitlines()}
            actual = [parse_rg_line(line) for line in input_rg.read_text().splitlines()]
            if len(set(actual)) != len(actual) or set(actual) != expected:
                raise ValueError(f"{input_rg}: @RG rows must match the observed IDs and their SM/LB/PU/PL/PM values")
            shutil.copyfile(input_rg, output_rg)
    return count


def parse_rg_line(line):
    kind, *fields = line.split("\t")
    if kind != "@RG" or len(fields) != 6 or any(":" not in tag for tag in fields):
        raise ValueError(f"invalid SAM read-group header row: {line!r}")
    return tuple(sorted(fields))


def write_manifest(spec, report_files, multiqc_files, read_groups):
    counts = {}
    for path in report_files:
        if path.name.endswith((".rna_read_retention.tsv", ".dna_read_retention.tsv")):
            with path.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if row["metric"] == "routed_branch_pairs":
                        key = (row["modality"], row["sample_id"], row["group"], row["branch"])
                        counts[key] = int(row["pairs"])

    def copy_metrics(paths, category):
        result = []
        for index, path in enumerate(sorted(paths)):
            dest = Path("preprocessing_qc") / category / str(index) / path.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, dest)
            result.append(str(Path("restart") / dest))
        return result

    headers = {path.name: path for path in read_groups}
    for split in spec["splits"]:
        key = (split["modality"], split["sample"], split["group"], split.get("mark", "RNA"))
        split["pairs"] = counts[key]
        if split["modality"] == "dna":
            published = Path(split["rg"])
            dest = Path("dna") / published.name
            dest.parent.mkdir(exist_ok=True)
            shutil.copyfile(headers[published.name], dest)
    spec["splits"].sort(key=lambda split: (split["modality"], split["split_id"]))
    report_paths = copy_metrics(report_files, "report")
    # Reuse identical report and MultiQC files in the published bundle.
    by_source = dict(zip(sorted(report_files), report_paths))
    extra = [path for path in multiqc_files if path not in by_source]
    spec["report_inputs"] = report_paths
    spec["multiqc_inputs"] = sorted([by_source[path] for path in multiqc_files if path in by_source]
                                    + copy_metrics(extra, "multiqc"))
    Path("restart_manifest.json").write_text(json.dumps(spec, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_subparsers(dest="mode", required=True)
    validate = modes.add_parser("validate")
    validate.add_argument("spec")
    validate.add_argument("r1", type=Path)
    validate.add_argument("r2", type=Path)
    validate.add_argument("--input-rg", type=Path)
    validate.add_argument("--output-rg", type=Path)
    write = modes.add_parser("write")
    write.add_argument("spec")
    write.add_argument("--report", nargs="*", type=Path, default=[])
    write.add_argument("--multiqc", nargs="*", type=Path, default=[])
    write.add_argument("--read-groups", nargs="*", type=Path, default=[])
    args = parser.parse_args()
    spec = json.loads(base64.b64decode(args.spec))
    try:
        if args.mode == "write":
            write_manifest(spec, args.report, args.multiqc, args.read_groups)
        else:
            count = validate_split(args.r1, args.r2, spec, args.input_rg, args.output_rg)
            Path("validated.txt").write_text(f"{count}\n")
    except (ValueError, OSError, KeyError) as error:
        split = spec.get("split_id", "manifest")
        if args.mode == "validate":
            split = f'{spec["modality"]}:{split} (R1: {spec.get("r1", args.r1)}, R2: {spec.get("r2", args.r2)})'
        print(f"ERROR: restart {split}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
