#!/usr/bin/env python3
"""Capture and compare semantic TrESFlow output contracts.

Binary genomics formats are decoded before comparison. Runtime-only reports
remain part of the path contract but their volatile contents are not compared.
"""

from __future__ import annotations

import argparse
import copy
import csv
import difflib
import gzip
import html
from html.parser import HTMLParser
import json
import math
import re
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any


CONTRACT_SCHEMA_VERSION = 1
ISO_TIMESTAMP = re.compile(
    r"\b\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:?\d{2})?\b"
)
STAR_TIMESTAMP = re.compile(r"\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\b")
ABSOLUTE_PATH = re.compile(r"(?<![A-Za-z0-9_.-])/(?:[^\s,;\"'<>()[\]{}]+/)*[^\s,;\"'<>()[\]{}]+")
CONTENT_PATH_ONLY = {
    "pipeline_info/execution_report.html",
    "pipeline_info/execution_timeline.html",
    "pipeline_info/execution_trace.tsv",
    "pipeline_info/flowchart.html",
}
STAR_MAPPING_SPEED_FIELD = "Mapping speed, Million of reads per hour"
STAR_MAPPING_SPEED_LINE = re.compile(
    rf"^(?P<prefix>\s*{re.escape(STAR_MAPPING_SPEED_FIELD)}\s*\|\s*)\S+\s*$"
)
PICARD_STARTED_ON_LINE = re.compile(r"^# Started on:\s+.+$")
RUNTIME_CONTRACT_PATH = "pipeline_info/runtime_contract.tsv"
RETIRED_PROCESS_HOST_TOOLS = {
    "codon",
    "cutadapt",
    "trim_galore",
    "pigz",
}


def normalize_scalar_text(text: str, roots: tuple[Path, ...] = ()) -> str:
    for root in sorted(roots, key=lambda item: len(str(item)), reverse=True):
        text = text.replace(str(root), "<PATH>")
    text = ISO_TIMESTAMP.sub("<TIMESTAMP>", text)
    text = STAR_TIMESTAMP.sub("<TIMESTAMP>", text)
    text = ABSOLUTE_PATH.sub("<PATH>", text)
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalize_lines(text: str, roots: tuple[Path, ...] = (), sort_lines: bool = False) -> list[str]:
    lines = [normalize_scalar_text(line.rstrip(), roots) for line in text.splitlines()]
    while lines and not lines[-1]:
        lines.pop()
    return sorted(lines) if sort_lines else lines


def normalize_scoped_runtime_line(relative_path: str, line: str) -> str:
    """Canonicalize only known runtime fields in their owning tool outputs."""
    if relative_path.endswith(".Log.final.out"):
        match = STAR_MAPPING_SPEED_LINE.fullmatch(line)
        if match:
            return f"{match.group('prefix')}<RUNTIME>"
    elif relative_path.endswith(".DuplicateMetrics.txt") and PICARD_STARTED_ON_LINE.fullmatch(line):
        return "# Started on: <TIMESTAMP>"
    return line


def canonicalize_contract_runtime_metadata(contract: dict[str, Any]) -> dict[str, Any]:
    """Return a comparison copy with narrowly scoped runtime metadata normalized."""
    canonical = copy.deepcopy(contract)
    for relative_path, lines in canonical.get("text", {}).items():
        canonical["text"][relative_path] = [
            normalize_scoped_runtime_line(relative_path, line) for line in lines
        ]
    runtime_contract = canonical.get("tables", {}).get(RUNTIME_CONTRACT_PATH)
    if runtime_contract:
        runtime_contract["rows"] = canonicalize_retired_process_host_rows(
            runtime_contract.get("rows", [])
        )
    return canonical


def canonicalize_retired_process_host_rows(rows: list[list[str]]) -> list[list[str]]:
    """Remove only host rows retired by process-environment migrations.

    Codon/Seq and FASTQ preprocessing now belong to process environments. The
    remaining host-tool rows and runtime paths are still compared as before.
    """
    canonical: list[list[str]] = []
    for row in rows:
        if row == ["[host_codon_seq_preflight]"]:
            break
        if row and row[0] in RETIRED_PROCESS_HOST_TOOLS | {"codon_home"}:
            continue
        canonical.append(row)
    return canonical


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def semantic_table(path: Path, roots: tuple[Path, ...]) -> dict[str, Any]:
    text = normalize_scalar_text(read_text(path), roots)
    lines = [line for line in text.splitlines() if line.strip()]
    delimiter = "\t" if any("\t" in line for line in lines[:2]) else ","
    rows = list(csv.reader(lines, delimiter=delimiter))
    return {"delimiter": "tab" if delimiter == "\t" else "comma", "rows": rows}


def semantic_matrix(path: Path) -> dict[str, Any]:
    header = ""
    dimensions: list[int] | None = None
    entries: list[list[str]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("%%MatrixMarket"):
                header = " ".join(line.split())
                continue
            if line.startswith("%"):
                continue
            fields = line.split()
            if dimensions is None:
                dimensions = [int(value) for value in fields]
            else:
                entries.append(fields)
    entries.sort(key=lambda row: tuple(int(value) if value.lstrip("-").isdigit() else value for value in row))
    return {"header": header, "dimensions": dimensions or [], "entries": entries}


def semantic_star_summary(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for row in csv.reader(read_text(path).splitlines()):
        if not row:
            continue
        rows.append([field.strip() for field in row])
    return rows


def canonical_header_line(line: str, roots: tuple[Path, ...]) -> str | None:
    fields = line.rstrip().split("\t")
    record_type = fields[0]
    tags = {}
    comments = []
    for field in fields[1:]:
        if ":" in field:
            key, value = field.split(":", 1)
            tags[key] = normalize_scalar_text(value, roots)
        else:
            comments.append(normalize_scalar_text(field, roots))
    if record_type == "@PG":
        # Program IDs, command lines, versions, and parent chains are generated
        # by tools and paths. Program identity is the stable semantic header.
        program = tags.get("PN") or tags.get("ID")
        return f"@PG\tPN:{program}" if program else None
    if record_type == "@HD":
        tags.pop("VN", None)
    if record_type in {"@HD", "@SQ", "@RG"}:
        return "\t".join([record_type, *[f"{key}:{tags[key]}" for key in sorted(tags)]])
    if record_type == "@CO":
        canonical = "\t".join(
            [record_type, *comments, *[f"{key}:{tags[key]}" for key in sorted(tags)]]
        )
        # STAR records its allocated thread count in @CO. Resource allocation
        # is runtime metadata and does not alter the alignment contract.
        return re.sub(r"(--runThreadN\s+)\d+", r"\1<THREADS>", canonical)
    return normalize_scalar_text(line.rstrip(), roots)


def semantic_bam(path: Path, samtools: Path, roots: tuple[Path, ...]) -> dict[str, Any]:
    header_result = subprocess.run(
        [str(samtools), "view", "-H", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    header = [
        canonical
        for line in header_result.stdout.splitlines()
        if (canonical := canonical_header_line(line, roots)) is not None
    ]
    # Dictionary order and @PG emission order are serialization details.
    header.sort(key=lambda line: (line.split("\t", 1)[0], line))

    view_result = subprocess.run(
        [str(samtools), "view", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    records: list[list[str]] = []
    for line in view_result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            raise ValueError(f"Malformed SAM record decoded from {path}: {line}")
        core = fields[:11]
        tags = sorted(fields[11:], key=lambda value: (value.split(":", 1)[0], value))
        records.append([*core, *tags])
    records.sort(key=lambda row: (row[0], int(row[1]), row[2], int(row[3]), row[5], row[6], int(row[7])))
    return {"header": header, "records": records}


def semantic_bigwig(path: Path) -> dict[str, Any]:
    try:
        import pyBigWig  # type: ignore[import-not-found]
    except ImportError as error:
        raise RuntimeError("pyBigWig is required for semantic bigWig regression checks") from error

    result: dict[str, Any] = {"chromosomes": {}, "intervals": {}}
    with pyBigWig.open(str(path)) as bigwig:
        chromosomes = bigwig.chroms()
        result["chromosomes"] = dict(sorted(chromosomes.items()))
        for chromosome in sorted(chromosomes):
            intervals = bigwig.intervals(chromosome)
            if not intervals:
                continue
            result["intervals"][chromosome] = [
                [int(start), int(end), round(float(value), 6)]
                for start, end, value in intervals
            ]
    return result


class VisibleHtmlParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.hidden_depth = 0
        self.text: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style"}:
            self.hidden_depth += 1

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style"} and self.hidden_depth:
            self.hidden_depth -= 1

    def handle_data(self, data: str) -> None:
        if not self.hidden_depth:
            value = " ".join(html.unescape(data).split())
            if value:
                self.text.append(value)


def semantic_html(path: Path, roots: tuple[Path, ...]) -> list[str]:
    parser = VisibleHtmlParser()
    parser.feed(read_text(path))
    visible = [normalize_scalar_text(value, roots) for value in parser.text]
    for index, value in enumerate(visible):
        if index and visible[index - 1] == "Pipeline version":
            visible[index] = "<PIPELINE_VERSION>"
    return visible


def normalize_json_value(value: Any, roots: tuple[Path, ...]) -> Any:
    if isinstance(value, dict):
        return {key: normalize_json_value(value[key], roots) for key in sorted(value)}
    if isinstance(value, list):
        return [normalize_json_value(item, roots) for item in value]
    if isinstance(value, str):
        return normalize_scalar_text(value, roots)
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return round(value, 10)
    return value


def semantic_json(path: Path, roots: tuple[Path, ...]) -> Any:
    return normalize_json_value(json.loads(read_text(path)), roots)


def semantic_multiqc_json(path: Path, roots: tuple[Path, ...]) -> dict[str, Any]:
    data = json.loads(read_text(path))
    stable_keys = (
        "report_general_stats_data",
        "report_general_stats_headers",
        "report_saved_raw_data",
    )
    return {key: normalize_json_value(data.get(key), roots) for key in stable_keys}


def semantic_fastqc_zip(path: Path, roots: tuple[Path, ...]) -> dict[str, list[str]]:
    contents: dict[str, list[str]] = {}
    with zipfile.ZipFile(path) as archive:
        for name in sorted(archive.namelist()):
            basename = Path(name).name
            if basename not in {"fastqc_data.txt", "summary.txt"}:
                continue
            text = archive.read(name).decode("utf-8", errors="replace")
            contents[basename] = normalize_lines(text, roots)
    if set(contents) != {"fastqc_data.txt", "summary.txt"}:
        raise ValueError(f"FastQC archive lacks semantic reports: {path}")
    return contents


def semantic_gzip_text(path: Path, roots: tuple[Path, ...]) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        # Tag-record order follows parallel serialization and is not biological.
        return normalize_lines(handle.read(), roots, sort_lines=True)


def capture_contract(root: Path, scenario: str, samtools: Path) -> dict[str, Any]:
    root = root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Output directory not found: {root}")
    roots = (root,)
    paths = sorted(path for path in root.rglob("*") if path.is_file())
    relative_files = [path.relative_to(root).as_posix() for path in paths]
    relative_directories = sorted(
        path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_dir()
    )

    contract: dict[str, Any] = {
        "scenario": scenario,
        "directories": relative_directories,
        "files": relative_files,
        "tables": {},
        "text": {},
        "json": {},
        "gzip_text": {},
        "star_matrices": {},
        "star_vectors": {},
        "star_summaries": {},
        "bams": {},
        "bigwigs": {},
        "html": {},
        "fastqc": {},
        "multiqc": {},
    }

    for path in paths:
        relative = path.relative_to(root).as_posix()
        if relative in CONTENT_PATH_ONLY:
            continue
        if path.suffix == ".bam":
            contract["bams"][relative] = semantic_bam(path, samtools, roots)
        elif path.suffix == ".bw":
            contract["bigwigs"][relative] = semantic_bigwig(path)
        elif path.name == "matrix.mtx" or path.name.endswith(".mtx"):
            contract["star_matrices"][relative] = semantic_matrix(path)
        elif path.name in {"barcodes.tsv", "features.tsv"} and "Solo.out" in relative:
            contract["star_vectors"][relative] = normalize_lines(read_text(path), roots)
        elif path.name == "Summary.csv" or path.name.endswith(".Solo.summary.csv"):
            contract["star_summaries"][relative] = semantic_star_summary(path)
        elif relative == "TrES_Stats/tres_report.html":
            contract["html"][relative] = semantic_html(path, roots)
        elif relative.endswith("multiqc_data.json"):
            contract["multiqc"][relative] = semantic_multiqc_json(path, roots)
        elif "/multiqc_report_data/" in relative and path.suffix == ".txt":
            contract["multiqc"][relative] = normalize_lines(read_text(path), roots)
        elif "/multiqc_report_data/" in relative:
            # Log, parquet, and citation serialization are runtime/package data.
            continue
        elif path.name.endswith("_fastqc.zip"):
            contract["fastqc"][relative] = semantic_fastqc_zip(path, roots)
        elif path.name.endswith("_fastqc.html") or relative.endswith("multiqc_report.html"):
            continue
        elif path.suffixes[-2:] == [".tsv", ".gz"]:
            contract["gzip_text"][relative] = semantic_gzip_text(path, roots)
        elif path.suffix == ".json":
            contract["json"][relative] = semantic_json(path, roots)
        elif path.suffix in {".tsv", ".csv"}:
            contract["tables"][relative] = semantic_table(path, roots)
        elif path.suffix in {".txt", ".stats", ".flagstat", ".idxstats", ".out", ".yml", ".yaml"}:
            contract["text"][relative] = normalize_lines(read_text(path), roots)
        elif path.suffix == ".bai" or path.suffix == ".parquet":
            continue
        else:
            # Genuinely deterministic small files only. Unknown large binary
            # formats must never silently fall back to byte comparison.
            if path.stat().st_size <= 4096:
                contract["text"][relative] = normalize_lines(read_text(path), roots)

    return canonicalize_contract_runtime_metadata(contract)


def validate_contract(contract: dict[str, Any]) -> None:
    """Reject incomplete captures before they can become accepted goldens."""
    scenario = contract["scenario"]
    errors: list[str] = []
    files = set(contract["files"])

    required_paths = {
        "TrES_Stats/barcode_composition.tsv",
        "TrES_Stats/library_complexity.tsv",
        "TrES_Stats/qc_metrics.tsv",
        "TrES_Stats/read_retention.tsv",
        "TrES_Stats/tres_report.html",
        "TrES_Stats/qc/multiqc/multiqc_report.html",
        "TrES_Stats/qc/multiqc/multiqc_report_data/multiqc_data.json",
        *CONTENT_PATH_ONLY,
    }
    missing_paths = sorted(required_paths - files)
    if missing_paths:
        errors.append("missing published paths: " + ", ".join(missing_paths))

    for section in ("tables", "gzip_text", "bams", "bigwigs", "html", "fastqc", "multiqc"):
        if not contract.get(section):
            errors.append(f"empty required semantic section: {section}")

    for path, bam in contract.get("bams", {}).items():
        if not bam["records"]:
            errors.append(f"BAM contains no alignment records: {path}")
    for path, bigwig in contract.get("bigwigs", {}).items():
        if not bigwig["intervals"]:
            errors.append(f"bigWig contains no intervals: {path}")

    if scenario == "rna_only":
        expected_counts = {
            "bams": 1,
            "bigwigs": 3,
            "fastqc": 3,
            "star_matrices": 3,
            "star_vectors": 4,
            "star_summaries": 1,
        }
        for section, expected_count in expected_counts.items():
            observed_count = len(contract.get(section, {}))
            if observed_count != expected_count:
                errors.append(
                    f"{section} count is {observed_count}; expected {expected_count} for {scenario}"
                )
        for path, matrix in contract.get("star_matrices", {}).items():
            if not matrix["dimensions"] or not matrix["entries"]:
                errors.append(f"STAR matrix is empty: {path}")
    elif scenario in {"dna_single", "dna_dual"}:
        expected_fastqc = 4 if scenario == "dna_single" else 3
        expected_counts = {"bams": 2, "bigwigs": 1, "fastqc": expected_fastqc}
        for section, expected_count in expected_counts.items():
            observed_count = len(contract.get(section, {}))
            if observed_count != expected_count:
                errors.append(
                    f"{section} count is {observed_count}; expected {expected_count} for {scenario}"
                )
        if not any(path.endswith(".DuplicateMetrics.txt") for path in contract.get("text", {})):
            errors.append("DNA duplicate metrics are missing")
        if not any(path.endswith(".dna_alignment_retention.tsv") for path in contract.get("tables", {})):
            errors.append("DNA alignment-retention metrics are missing")
        artifact_paths = [
            path
            for path in contract.get("json", {})
            if path.endswith(".dual_tag_artifact_filter.cutadapt.json")
        ]
        if scenario == "dna_dual" and not artifact_paths:
            errors.append("dual-tag artifact-filter report is missing")
        elif scenario == "dna_dual":
            report = contract["json"][artifact_paths[0]]
            read_counts = report.get("read_counts", {})
            filtered = read_counts.get("filtered", {})
            observed = {
                "input": read_counts.get("input"),
                "output": read_counts.get("output"),
                "discard_trimmed": filtered.get("discard_trimmed"),
            }
            expected = {"input": 64, "output": 63, "discard_trimmed": 1}
            if observed != expected:
                errors.append(
                    f"dual-tag artifact-filter counts are {observed}; expected {expected}"
                )
    else:
        errors.append(f"unknown scenario: {scenario}")

    if errors:
        raise ValueError("Incomplete semantic regression contract:\n- " + "\n- ".join(errors))


def remove_empty_sections(contract: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in contract.items() if value not in ({}, [])}


def write_capture(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite an existing baseline contract: {output}")
    provenance = json.loads(args.provenance.read_text(encoding="utf-8"))
    contract = capture_contract(args.root, args.scenario, args.samtools)
    validate_contract(contract)
    contract = remove_empty_sections(contract)
    payload = {
        "contract_schema_version": CONTRACT_SCHEMA_VERSION,
        "baseline_provenance": provenance,
        "contract": contract,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def compare_capture(args: argparse.Namespace) -> None:
    expected_payload = json.loads(args.expected.read_text(encoding="utf-8"))
    expected = canonicalize_contract_runtime_metadata(expected_payload["contract"])
    validate_contract(expected)
    actual = canonicalize_contract_runtime_metadata(
        capture_contract(args.root, expected["scenario"], args.samtools)
    )
    validate_contract(actual)
    actual = remove_empty_sections(actual)
    if actual == expected:
        print(f"PASS: semantic contract matches {args.expected}")
        return
    expected_text = json.dumps(expected, indent=2, sort_keys=True).splitlines()
    actual_text = json.dumps(actual, indent=2, sort_keys=True).splitlines()
    diff = "\n".join(
        difflib.unified_diff(expected_text, actual_text, fromfile=str(args.expected), tofile=str(args.root), n=3)
    )
    print(diff[:100_000], file=sys.stderr)
    raise SystemExit("Semantic regression contract mismatch")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samtools", type=Path, default=Path("samtools"), help="samtools executable used to decode BAM"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    capture = subparsers.add_parser("capture")
    capture.add_argument("--root", required=True, type=Path)
    capture.add_argument("--scenario", required=True, choices=("rna_only", "dna_single", "dna_dual"))
    capture.add_argument("--provenance", required=True, type=Path)
    capture.add_argument("--output", required=True, type=Path)
    capture.set_defaults(handler=write_capture)

    compare = subparsers.add_parser("compare")
    compare.add_argument("--root", required=True, type=Path)
    compare.add_argument("--expected", required=True, type=Path)
    compare.set_defaults(handler=compare_capture)
    return parser.parse_args()


if __name__ == "__main__":
    parsed_args = parse_args()
    parsed_args.handler(parsed_args)
