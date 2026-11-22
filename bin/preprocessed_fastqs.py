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

from run_fq_to_sam import extract_rna_tags, normalize_qname, parse_header
from tresflow_fastq_utils import fastq_iter, parse_aviti_qname, write_rg_header


BULK_RNA_SUFFIX = re.compile(r"_P[0-9a-f]{64}_[^_]+$")
LIGATION_BARCODE = re.compile(r"[ACGT]+")
UMI_BASES = re.compile(r"[ACGTN]+")
UMI_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
REQUIRED_TAGS = {"CB", "XI", "SB", "RG", "UM", "MO"}


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
    if not dna and BULK_RNA_SUFFIX.search(name):
        raise ValueError("RNA bulk export detected. Supply the original, unclipped post-split RNA FASTQs.")

    tags = read_tags(comment1)
    required = {"CB", "XI", "SB", "RG", "MO" if dna else "UM"}
    if not required <= tags.keys():
        raise ValueError(f"missing required FASTQ tags: {', '.join(sorted(required - tags.keys()))}")
    if tags != read_tags(comment2):
        raise ValueError("mate FASTQ tags differ")
    sb_index = barcodes.get(tags["SB"])
    prefix = f'{spec["sample"]}_{spec["group"]}_{sb_index}_'
    cell = tags["CB"]
    ligation = cell.removeprefix(prefix)
    if (sb_index is None or not cell.startswith(prefix) or cell != tags["XI"]
            or len(ligation) != ligation_length or not LIGATION_BARCODE.fullmatch(ligation)):
        raise ValueError("CB/XI and corrected SB do not match the samplesheet group and logical barcode")

    if dna:
        # The Codon splitter preserves the lane token, including any leading zeros.
        parse_aviti_qname(name)
        instrument, run, flowcell, lane, *_ = name.split(":")
        expected_rg = f"{instrument}:{run}:{flowcell}:L{lane}"
        if tags["RG"] != expected_rg or tags["MO"] != spec["mo_barcode"]:
            raise ValueError("RG does not match AVITI QNAME or MO does not match the DNA mark")
        return expected_rg

    if tags["RG"] != cell or len(tags["UM"]) != spec["umi_length"]:
        raise ValueError("RNA RG must equal CB/XI and UM must match the configured UMI length")
    umi_start = spec["umi_start"]
    umi_end = umi_start + spec["umi_length"]
    if len(second[1]) >= umi_end:
        bases = second[1][umi_start:umi_end].upper()
        if UMI_BASES.fullmatch(bases) and bases.translate(UMI_COMPLEMENT)[::-1] != tags["UM"].upper():
            raise ValueError("RNA R2 prefix does not encode the UM tag. Supply unclipped post-split FASTQs.")
    if extract_rna_tags(comment1) != extract_rna_tags(comment2):
        raise ValueError("mate RNA SAM tags differ")
    return None


def validate_split(r1, r2, spec, input_rg=None, output_rg=None):
    count = 0
    groups = set()
    ligation_length = spec["ligation_length"]
    barcodes = spec["sb_indices"]
    try:
        for count, (first, second) in enumerate(zip_longest(fastq_iter(r1), fastq_iter(r2)), 1):
            if first is None or second is None:
                raise ValueError("mate FASTQs have unequal record counts")
            rg = validate_pair(first, second, spec, barcodes, ligation_length)
            if rg is not None:
                groups.add(rg)
    except (ValueError, OSError, EOFError) as error:
        raise ValueError(f"at pair {count or 1}: {error}") from error

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
        print(f'ERROR: restart {spec.get("split_id", "manifest")}: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
