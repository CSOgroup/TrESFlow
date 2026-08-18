"""Pure parsers, validation, and normalized report model for TrESFlow QC."""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional


RETENTION_FIELDS = [
    "run", "parent_sample", "modality", "branch", "stage_order", "stage",
    "stage_scope", "input_pairs", "output_pairs", "retained_prev_pct",
    "cumulative_raw_pct", "unit", "count_source", "subset_verified",
]
QC_FIELDS = [
    "run", "parent_sample", "modality", "branch", "metric", "value",
    "unit", "denominator", "source",
]
COMPOSITION_FIELDS = [
    "run", "parent_sample", "modality", "group", "barcode_type", "category",
    "barcode_sequence", "barcode_label", "count", "percentage",
    "denominator_count", "denominator_definition", "contract_level", "source",
]
COMPLEXITY_FIELDS = [
    "run", "parent_sample", "branch", "read_pairs_examined",
    "observed_unique_pairs", "pcr_or_library_duplicate_pairs",
    "optical_duplicate_pairs", "observed_unique_fraction",
    "pcr_or_library_duplicate_fraction", "optical_duplicate_fraction",
    "percent_duplication", "estimated_library_size", "roi_depth_multiplier",
    "roi_yield_multiplier", "roi_estimated_unique_pairs", "source",
]


@dataclass(frozen=True)
class Branch:
    run: str
    parent: str
    modality: str
    branch: str
    group: str
    split_id: str


@dataclass
class ReportModel:
    library_name: str
    pipeline_version: str
    generated_at: str
    branches: list[Branch]
    retention: list[dict[str, object]]
    qc_metrics: list[dict[str, object]]
    barcode_composition: list[dict[str, object]]
    library_complexity: list[dict[str, object]]
    insights: list[dict[str, object]]
    warnings: list[str]
    run_metadata: dict[str, object] = field(default_factory=dict)
    raw: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": "2.0",
            "library_name": self.library_name,
            "pipeline_version": self.pipeline_version,
            "generated_at": self.generated_at,
            "branches": [asdict(branch) for branch in self.branches],
            "retention": self.retention,
            "qc_metrics": self.qc_metrics,
            "barcode_composition": self.barcode_composition,
            "library_complexity": self.library_complexity,
            "insights": self.insights,
            "warnings": self.warnings,
            "run_metadata": self.run_metadata,
            # Compatibility for existing consumers of tres_report_metrics.json.
            "raw": self.raw,
        }


def retention_display_stage(stage: object) -> str:
    """Return the report-facing label without changing the accounting key."""
    value = str(stage)
    return {
        "Properly paired": "Mapped",
        "Proper pair": "Mapped",
        "Properly paired after blacklist": "Mapped",
    }.get(value, value)


def numeric(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text.endswith("%"):
        text = text[:-1]
    if text.lower() in {"", "nan", "-nan", "n/a", "na", "none"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def integer(value: object) -> Optional[int]:
    parsed = numeric(value)
    return int(round(parsed)) if parsed is not None else None


def percentage(numerator: Optional[int], denominator: Optional[int]) -> Optional[float]:
    if numerator is None or denominator in {None, 0}:
        return None
    return numerator / denominator * 100.0


def fraction(numerator: Optional[int], denominator: Optional[int]) -> Optional[float]:
    if numerator is None or denominator in {None, 0}:
        return None
    return numerator / denominator


def find_files(root: Path) -> list[Path]:
    """Find report inputs without walking BAM/FASTQ-heavy result trees."""
    if root.name == "inputs":
        return [path for path in root.rglob("*") if path.is_file()]
    patterns = (
        "TrES_Stats/*.tsv",
        "TrES_Stats/qc/samtools/*.flagstat",
        "dna_align/*.DuplicateMetrics.txt",
        "rna_align/*.Log.final.out",
        "rna_align/*.Solo.outGeneFull/Summary.csv",
        "pipeline_info/derived_contract/*.tsv",
    )
    found: list[Path] = []
    for pattern in patterns:
        found.extend(path for path in root.glob(pattern) if path.is_file())
    # A flat or synthetic fixture remains a supported standalone input form.
    if not found:
        found.extend(path for path in root.rglob("*") if path.is_file())
    return list(dict.fromkeys(found))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_stats(path: Path) -> dict[str, dict[str, Optional[float]]]:
    result: dict[str, dict[str, Optional[float]]] = {}
    with path.open(errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if not fields or not fields[0] or fields[0].startswith("#"):
                continue
            result[fields[0]] = {
                "count": integer(fields[1]) if len(fields) > 1 else None,
                "percent": numeric(fields[2]) if len(fields) > 2 else None,
            }
    return result


def read_barcode_counts(path: Path) -> dict[str, int]:
    """Read the headerless ``count<TAB>barcode`` contract emitted by Tag.

    A headered ``barcode<TAB>count`` form is accepted as well so fixtures and
    future wrappers can make the schema explicit without breaking the reader.
    Repeated barcode rows are summed; negative or malformed counts fail rather
    than silently changing the composition denominator.
    """
    result: dict[str, int] = defaultdict(int)
    with path.open(errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0] or fields[0].startswith("#"):
                continue
            count = integer(fields[0])
            barcode = fields[1].strip()
            if count is None:
                # Optional explicit header or barcode-first compatibility.
                reverse_count = integer(fields[1])
                if fields[0].strip().lower() in {"barcode", "barcode_sequence", "sequence"}:
                    continue
                if reverse_count is not None:
                    barcode, count = fields[0].strip(), reverse_count
            if count is None or count < 0 or not barcode:
                raise ValueError(f"Malformed sample-barcode count at {path}:{line_number}")
            result[barcode] += int(count)
    return dict(result)


def read_star_log(path: Path) -> dict[str, Optional[float]]:
    result: dict[str, Optional[float]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if "|" in line:
                key, value = line.split("|", 1)
                result[key.strip()] = numeric(value)
    return result


def read_summary(path: Path) -> dict[str, Optional[float]]:
    result: dict[str, Optional[float]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if "," in line:
                key, value = line.rstrip("\n").split(",", 1)
                result[key.strip()] = numeric(value)
    return result


def read_flagstat(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    patterns = {
        "total": r"^(\d+) \+ \d+ in total",
        "primary": r"^(\d+) \+ \d+ primary$",
        "mapped": r"^(\d+) \+ \d+ mapped ",
        "primary_mapped": r"^(\d+) \+ \d+ primary mapped ",
        "paired": r"^(\d+) \+ \d+ paired in sequencing",
        "read1": r"^(\d+) \+ \d+ read1",
        "read2": r"^(\d+) \+ \d+ read2",
        "properly_paired": r"^(\d+) \+ \d+ properly paired ",
        "duplicates": r"^(\d+) \+ \d+ duplicates",
    }
    with path.open(errors="replace") as handle:
        for line in handle:
            for name, pattern in patterns.items():
                match = re.match(pattern, line)
                if match:
                    result[name] = int(match.group(1))
                    break
    return result


def read_duplicate_metrics(path: Path) -> dict[str, object]:
    """Parse Picard's metrics row and actual ROI histogram points."""
    header: Optional[list[str]] = None
    histogram_header: Optional[list[str]] = None
    metrics: dict[str, object] = {"file": path.name}
    roi: list[dict[str, float]] = []
    with path.open(errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                header = None if len(metrics) > 1 else header
                histogram_header = None
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if "PERCENT_DUPLICATION" in fields:
                header = fields
                continue
            if header is not None and "read_pairs_examined" not in metrics:
                record = dict(zip(header, fields))
                duplication_fraction = numeric(record.get("PERCENT_DUPLICATION"))
                metrics.update(
                    {
                        "library": record.get("LIBRARY", ""),
                        "read_pairs_examined": integer(record.get("READ_PAIRS_EXAMINED")),
                        "read_pair_duplicates": integer(record.get("READ_PAIR_DUPLICATES")),
                        "read_pair_optical_duplicates": integer(record.get("READ_PAIR_OPTICAL_DUPLICATES")),
                        "percent_duplication": duplication_fraction * 100.0 if duplication_fraction is not None else None,
                        "estimated_library_size": integer(record.get("ESTIMATED_LIBRARY_SIZE")),
                    }
                )
                header = None
                continue
            if "BIN" in fields and "CoverageMult" in fields:
                histogram_header = fields
                continue
            if histogram_header is not None:
                record = dict(zip(histogram_header, fields))
                depth = numeric(record.get("BIN"))
                yield_multiplier = numeric(record.get("CoverageMult"))
                # The same Picard file can contain duplicate-set bins beyond
                # the ROI domain. CoverageMult ROI is defined for 1..100.
                if depth is not None and 1 <= depth <= 100 and yield_multiplier is not None:
                    roi.append({"depth_multiplier": float(depth), "yield_multiplier": float(yield_multiplier)})
    if "read_pairs_examined" in metrics:
        metrics["roi"] = roi
    return metrics


def read_artifact_summary(path: Path) -> dict[str, object]:
    rows = read_tsv(path)
    if len(rows) != 1:
        raise ValueError(f"Expected one dual-tag artifact summary row in {path}")
    row = rows[0]
    names = ("input_pairs", "retained_pairs", "rejected_pairs", "r1_with_signature", "r2_with_signature")
    values = {name: integer(row.get(name)) for name in names}
    if any(value is None or value < 0 for value in values.values()):
        raise ValueError(f"Invalid non-negative integer in dual-tag artifact summary {path}")
    input_pairs = int(values["input_pairs"])
    retained = int(values["retained_pairs"])
    rejected = int(values["rejected_pairs"])
    r1 = int(values["r1_with_signature"])
    r2 = int(values["r2_with_signature"])
    if input_pairs != retained + rejected:
        raise ValueError(f"Dual-tag artifact accounting mismatch in {path}: {input_pairs} != {retained} + {rejected}")
    both = r1 + r2 - rejected
    r1_only = r1 - both
    r2_only = r2 - both
    if min(both, r1_only, r2_only) < 0:
        raise ValueError(f"Invalid mate-overlap arithmetic in dual-tag artifact summary {path}")
    return {
        **values,
        "both_mates": both,
        "r1_only": r1_only,
        "r2_only": r2_only,
        "cutadapt_version": row.get("cutadapt_version", ""),
        "signature_fasta_sha256": row.get("signature_fasta_sha256", ""),
        "sample_id": row.get("sample_id", ""),
        "tagmentation": row.get("tagmentation", "dual"),
        "file": path.name,
        "path": path.name,
    }


def _contract_rows(files: list[Path], name: str) -> list[dict[str, str]]:
    candidates = [path for path in files if path.name == name]
    if not candidates:
        return []
    candidates.sort(key=lambda path: ("derived_contract" not in path.parts, len(path.parts)))
    return read_tsv(candidates[0])


def _unique_run_labels(pairs: Iterable[tuple[str, str]]) -> dict[tuple[str, str], str]:
    pairs = list(dict.fromkeys(pairs))
    owners: dict[str, set[str]] = defaultdict(set)
    for parent, group in pairs:
        owners[group].add(parent)
    return {(parent, group): group if len(owners[group]) == 1 else f"{parent}_{group}" for parent, group in pairs}


def discover_branches(files: list[Path]) -> list[Branch]:
    rna_pairs: list[tuple[str, str]] = []
    dna_specs: list[tuple[str, str, str]] = []
    group_pairs: list[tuple[str, str]] = []

    for path in files:
        match = re.match(r"^(.+)\.rna_read_retention\.tsv$", path.name)
        if match:
            parent = match.group(1)
            for row in read_tsv(path):
                if row.get("metric") == "routed_branch_pairs":
                    pair = (parent, row.get("group", ""))
                    if all(pair) and pair not in rna_pairs:
                        rna_pairs.append(pair)
                        group_pairs.append(pair)
        match = re.match(r"^(.+)\.dna_read_retention\.tsv$", path.name)
        if match:
            parent = match.group(1)
            for row in read_tsv(path):
                if row.get("metric") == "routed_branch_pairs":
                    spec = (parent, row.get("group", ""), row.get("branch", ""))
                    if all(spec) and spec not in dna_specs:
                        dna_specs.append(spec)
                        group_pairs.append(spec[:2])

    for row in _contract_rows(files, "rna_sb_group_map.tsv"):
        pair = (row.get("sample", ""), row.get("sb_group", ""))
        if all(pair) and pair not in rna_pairs:
            rna_pairs.append(pair)
            group_pairs.append(pair)
    for row in _contract_rows(files, "dna_mo_map.tsv"):
        spec = (row.get("sample", ""), row.get("sb_group", ""), row.get("mark", ""))
        if all(spec) and spec not in dna_specs:
            dna_specs.append(spec)
            group_pairs.append(spec[:2])

    # Last-resort legacy discovery. Exact group identity cannot be recovered
    # from filename tokens in every naming scheme, so mark this as coarse.
    if not rna_pairs:
        for path in files:
            if path.name == "Summary.csv" and path.parent.name.endswith(".Solo.outGeneFull"):
                split_id = path.parent.name.removesuffix(".Solo.outGeneFull")
                rna_pairs.append((split_id, split_id))
                group_pairs.append((split_id, split_id))
    if not dna_specs:
        for path in files:
            match = re.match(r"^dna\.(.+)\.nodup\.flagstat$", path.name)
            if match:
                split_id = match.group(1)
                parent_group, _, mark = split_id.rpartition("_")
                spec = (parent_group or split_id, parent_group or split_id, mark or "DNA")
                if spec not in dna_specs:
                    dna_specs.append(spec)
                    group_pairs.append(spec[:2])

    labels = _unique_run_labels(group_pairs)
    branches: list[Branch] = []
    for parent, group in rna_pairs:
        split_id = parent if parent == group else f"{parent}_{group}"
        branches.append(Branch(labels[(parent, group)], parent, "RNA", "RNA", group, split_id))
    for parent, group, mark in dna_specs:
        split_id = f"{parent}_{mark}" if parent == group else f"{parent}_{group}_{mark}"
        branches.append(Branch(labels[(parent, group)], parent, "DNA", mark, group, split_id))
    return sorted(branches, key=lambda value: (value.run, value.modality != "RNA", value.branch))


def _retention_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], int]:
    result = {}
    for row in rows:
        count = integer(row.get("pairs"))
        if count is not None:
            result[(row.get("group", ""), row.get("branch", ""), row.get("metric", ""))] = count
    return result


def _metric_lookup(rows: list[dict[str, str]]) -> dict[str, int]:
    return {row.get("metric", ""): int(integer(row.get("pairs"))) for row in rows if row.get("metric") and integer(row.get("pairs")) is not None}


def _gate_lookup(rows: list[dict[str, str]], modality: str) -> dict[str, int]:
    values = _metric_lookup(rows)
    expected = ["split_input_pairs", "ligation_barcode_accepted_pairs", "sample_barcode_accepted_pairs"]
    if modality == "dna":
        expected.append("modality_barcode_accepted_pairs")
    missing = [name for name in expected if name not in values]
    if missing:
        raise ValueError(f"Incomplete exact {modality.upper()} barcode-gate contract: {missing}")
    counts = [values[name] for name in expected]
    if any(right > left for left, right in zip(counts, counts[1:])):
        raise ValueError(f"Exact {modality.upper()} barcode gates are not nested: {values}")
    return values


def _raw_pairs(stats: dict[tuple[str, str], tuple[dict, Path]], parent: str, modality: str) -> tuple[Optional[int], str]:
    entry = stats.get((parent, modality.lower()))
    if not entry:
        return None, ""
    values, path = entry
    reads = values.get("reads", {}).get("count")
    return integer(reads), path.name


def _add_path(destination: list[dict[str, object]], branch: Branch, raw: int, raw_source: str, stages: list[tuple]) -> None:
    previous = raw
    all_stages = [("Raw input", raw, "shared", raw_source, "read pairs")] + stages
    for order, (name, count, scope, source, unit) in enumerate(all_stages, start=1):
        count = int(count)
        if count < 0 or count > previous:
            raise ValueError(f"Non-nested retention path for {branch.split_id}: {name}={count:,}, preceding={previous:,}")
        destination.append(
            {
                "run": branch.run, "parent_sample": branch.parent, "modality": branch.modality,
                "branch": branch.branch, "stage_order": order, "stage": name,
                "stage_scope": scope, "input_pairs": previous, "output_pairs": count,
                "retained_prev_pct": percentage(count, previous), "cumulative_raw_pct": percentage(count, raw),
                "unit": unit, "count_source": source, "subset_verified": "yes",
            }
        )
        previous = count


def _add_qc(rows: list[dict[str, object]], branch: Branch, metric: str, value: object, unit: str, denominator: str, source: str, branch_name: Optional[str] = None) -> None:
    if value is not None:
        rows.append(
            {
                "run": branch.run, "parent_sample": branch.parent, "modality": branch.modality,
                "branch": branch.branch if branch_name is None else branch_name,
                "metric": metric, "value": value, "unit": unit,
                "denominator": denominator, "source": source,
            }
        )


def _parse_composition(files: list[Path], branches: list[Branch], warnings: list[str]) -> list[dict[str, object]]:
    run_for_group = {(branch.parent, branch.group): branch.run for branch in branches}
    result: list[dict[str, object]] = []
    legacy_sample_contracts: list[str] = []
    for path in files:
        match = re.match(r"^(.+)\.(rna|dna)_barcode_composition\.tsv$", path.name)
        if not match:
            continue
        parent, modality = match.groups()
        rows = read_tsv(path)
        sample_rows = [row for row in rows if row.get("barcode_type") == "sample_barcode"]
        mark_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            if row.get("barcode_type") == "dna_mark":
                mark_groups[row.get("group", "__all__")].append(row)
        category_sets = [("sample_barcode", "__all__", sample_rows)] + [("dna_mark", group, values) for group, values in mark_groups.items()]
        for barcode_type, logical_group, category_rows in category_sets:
            if not category_rows:
                continue
            denominators = {integer(row.get("denominator_count")) for row in category_rows}
            if None in denominators or len(denominators) != 1:
                raise ValueError(f"Inconsistent composition denominator in {path}: {barcode_type}/{logical_group}")
            denominator = int(next(iter(denominators)))
            counts = [integer(row.get("count")) for row in category_rows]
            if any(count is None or count < 0 for count in counts) or sum(int(count) for count in counts) != denominator:
                raise ValueError(f"Composition categories do not reconcile in {path}: {barcode_type}/{logical_group}")
            if not any(row.get("category") == "NoMatch" for row in category_rows):
                raise ValueError(f"Composition contract omits NoMatch in {path}: {barcode_type}/{logical_group}")
            is_sequence_level = barcode_type != "sample_barcode" or all(
                row.get("category") == "NoMatch" or bool(row.get("barcode_sequence")) for row in category_rows
            )
            contract_level = "exact_sequence" if is_sequence_level else "legacy_group"
            if barcode_type == "sample_barcode" and not is_sequence_level:
                legacy_sample_contracts.append(f"{parent} {modality.upper()}")
            parent_runs = {branch.run for branch in branches if branch.parent == parent}
            sample_composition_run = next(iter(parent_runs)) if len(parent_runs) == 1 else parent
            for row in category_rows:
                category = row.get("category", "")
                group = row.get("group", logical_group)
                if barcode_type == "sample_barcode" and group == "__all__" and category != "NoMatch":
                    group = category if not is_sequence_level else group
                run = sample_composition_run if barcode_type == "sample_barcode" else run_for_group.get((parent, group))
                if not run:
                    run = sample_composition_run
                count = int(integer(row.get("count")))
                result.append(
                    {
                        "run": run, "parent_sample": parent, "modality": modality.upper(),
                        "group": group, "barcode_type": barcode_type, "category": category,
                        "barcode_sequence": row.get("barcode_sequence", ""),
                        "barcode_label": row.get("barcode_label", "") or category,
                        "count": count, "percentage": percentage(count, denominator) if denominator else 0.0,
                        "denominator_count": denominator,
                        "denominator_definition": row.get("denominator_definition", ""),
                        "contract_level": contract_level, "source": path.name,
                    }
                )
    return result


def _sample_barcode_contracts(files: list[Path], modality: str) -> dict[tuple[str, str], dict[str, str]]:
    """Return configured sequence metadata keyed by ``(sample, sequence)``."""
    rows = _contract_rows(files, f"{modality}_sb_group_map.tsv")
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        parent = row.get("sample", "").strip()
        sequence = row.get("sb_bc", "").strip()
        group = row.get("sb_group", "").strip()
        if not parent or not sequence or not group:
            continue
        key = (parent, sequence)
        previous = result.get(key)
        if previous and previous["group"] != group:
            raise ValueError(
                f"Conflicting {modality.upper()} sample-barcode groups for {parent}/{sequence}: "
                f"{previous['group']} vs {group}"
            )
        label = (
            row.get("barcode_label", "").strip()
            or row.get("sb_label", "").strip()
            or row.get("label", "").strip()
            or sequence
        )
        result[key] = {"group": group, "label": label}
    return result


def _sample_composition_from_counts(
    files: list[Path],
    branches: list[Branch],
) -> tuple[list[dict[str, object]], set[tuple[str, str]]]:
    """Build exact per-sequence sample composition from Tag counts/stats.

    Assigned sequence counts come only from ``*.counts.tsv``.  ``NoMatch`` is
    taken from the paired stats file's ``reads_without_bc`` metric.  The reader
    validates both ``bc_reads`` and ``reads`` when present, including any
    redundant NoMatch row retained in historical counts files.
    """
    paths = {path.name: path for path in files}
    run_for_group = {(branch.parent, branch.group): branch.run for branch in branches}
    parent_runs: dict[str, set[str]] = defaultdict(set)
    for branch in branches:
        parent_runs[branch.parent].add(branch.run)
    rows: list[dict[str, object]] = []
    covered: set[tuple[str, str]] = set()

    for counts_path in files:
        match = re.match(r"^(.+)\.(rna|dna)_sample_barcode\.counts\.tsv$", counts_path.name)
        if not match:
            continue
        parent, modality = match.groups()
        stats_path = paths.get(f"{parent}.{modality}_sample_barcode.stats.tsv")
        if stats_path is None:
            raise ValueError(f"Sample-barcode counts lack matching stats file: {counts_path.name}")
        stats = read_stats(stats_path)
        counts = read_barcode_counts(counts_path)
        contracts = _sample_barcode_contracts(files, modality)
        configured = {
            sequence: metadata
            for (sample, sequence), metadata in contracts.items()
            if sample == parent
        }
        assigned_counts = {
            sequence: int(count)
            for sequence, count in counts.items()
            if sequence.lower() != "nomatch"
        }
        if configured:
            unexpected = sorted(set(assigned_counts) - set(configured))
            if unexpected:
                raise ValueError(
                    f"Unconfigured {modality.upper()} sample barcode(s) in {counts_path.name}: "
                    + ", ".join(unexpected)
                )
            for sequence in configured:
                assigned_counts.setdefault(sequence, 0)
        assigned = sum(assigned_counts.values())
        stats_assigned = integer(stats.get("bc_reads", {}).get("count"))
        no_match = integer(stats.get("reads_without_bc", {}).get("count"))
        total = integer(stats.get("reads", {}).get("count"))
        if no_match is None:
            raise ValueError(f"Missing reads_without_bc in {stats_path.name}")
        redundant_no_match = next(
            (int(value) for key, value in counts.items() if key.lower() == "nomatch"),
            None,
        )
        if redundant_no_match is not None and redundant_no_match != no_match:
            raise ValueError(
                f"NoMatch disagreement for {parent} {modality.upper()}: "
                f"{counts_path.name}={redundant_no_match} vs {stats_path.name}={no_match}"
            )
        if stats_assigned is not None and assigned != stats_assigned:
            raise ValueError(
                f"Assigned sample-barcode counts disagree for {parent} {modality.upper()}: "
                f"{assigned} != bc_reads {stats_assigned}"
            )
        denominator = assigned + no_match
        if total is not None and denominator != total:
            raise ValueError(
                f"Sample-barcode composition does not reconcile for {parent} {modality.upper()}: "
                f"assigned {assigned} + reads_without_bc {no_match} != reads {total}"
            )
        definition = (
            f"{modality.upper()} pairs evaluated by sample-barcode tagging "
            "(configured-barcode assignments plus reads_without_bc)"
        )
        unique_parent_run = next(iter(parent_runs[parent])) if len(parent_runs[parent]) == 1 else parent
        for sequence, count in sorted(assigned_counts.items()):
            metadata = configured.get(sequence, {"group": "__all__", "label": sequence})
            group = metadata["group"]
            rows.append(
                {
                    "run": run_for_group.get((parent, group), unique_parent_run),
                    "parent_sample": parent,
                    "modality": modality.upper(),
                    "group": group,
                    "barcode_type": "sample_barcode",
                    "category": sequence,
                    "barcode_sequence": sequence,
                    "barcode_label": metadata["label"],
                    "count": count,
                    "percentage": percentage(count, denominator) if denominator else 0.0,
                    "denominator_count": denominator,
                    "denominator_definition": definition,
                    "contract_level": "exact_sequence_counts_stats",
                    "source": f"{counts_path.name};{stats_path.name}",
                }
            )
        rows.append(
            {
                "run": unique_parent_run,
                "parent_sample": parent,
                "modality": modality.upper(),
                "group": "__all__",
                "barcode_type": "sample_barcode",
                "category": "NoMatch",
                "barcode_sequence": "",
                "barcode_label": "NoMatch",
                "count": no_match,
                "percentage": percentage(no_match, denominator) if denominator else 0.0,
                "denominator_count": denominator,
                "denominator_definition": definition,
                "contract_level": "exact_sequence_counts_stats",
                "source": f"{counts_path.name};{stats_path.name}:reads_without_bc",
            }
        )
        covered.add((parent, modality))
    return rows, covered


def _collect_marginal_qc(qc: list[dict[str, object]], branch: Branch, values: dict, path: Path, prefix: str) -> None:
    matched = values.get("bc_reads", {})
    _add_qc(qc, branch, f"{prefix}_matched_pairs", matched.get("count"), "read_pairs", "raw input pairs; marginal", path.name, "ALL")
    _add_qc(qc, branch, f"{prefix}_match_pct", matched.get("percent"), "percent", "raw input pairs; marginal", path.name, "ALL")


def _validate_composition_against_gates(
    composition: list[dict[str, object]],
    gates: dict[tuple[str, str], tuple[list[dict[str, str]], Path]],
) -> None:
    """Reconcile exhaustive categories to the exact cumulative gate contract."""
    by_parent_modality: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in composition:
        by_parent_modality[(str(row["parent_sample"]), str(row["modality"]).lower())].append(row)
    for key, rows in by_parent_modality.items():
        gate_contract = gates.get(key)
        if not gate_contract:
            continue
        values = _gate_lookup(gate_contract[0], key[1])
        sample_rows = [row for row in rows if row["barcode_type"] == "sample_barcode"]
        if sample_rows:
            denominators = {int(row["denominator_count"]) for row in sample_rows}
            accepted = sum(int(row["count"]) for row in sample_rows if row["category"] != "NoMatch")
            if denominators != {values["ligation_barcode_accepted_pairs"]}:
                raise ValueError(f"Sample-barcode composition denominator disagrees with exact ligation gate for {key}")
            if accepted != values["sample_barcode_accepted_pairs"]:
                raise ValueError(f"Sample-barcode categories disagree with exact accepted gate for {key}: {accepted} != {values['sample_barcode_accepted_pairs']}")
        if key[1] != "dna":
            continue
        mark_rows = [row for row in rows if row["barcode_type"] == "dna_mark"]
        by_group: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in mark_rows:
            by_group[str(row["group"])].append(row)
        sample_by_group: dict[str, int] = defaultdict(int)
        for row in sample_rows:
            if row["category"] != "NoMatch":
                sample_by_group[str(row["group"])] += int(row["count"])
        for group, group_rows in by_group.items():
            denominators = {int(row["denominator_count"]) for row in group_rows}
            if denominators != {sample_by_group[group]}:
                raise ValueError(f"DNA mark-composition denominator disagrees with sample assignments for {key}/{group}")
        mark_accepted = sum(int(row["count"]) for row in mark_rows if row["category"] != "NoMatch")
        if mark_accepted != values["modality_barcode_accepted_pairs"]:
            raise ValueError(f"DNA mark categories disagree with exact MO gate for {key}: {mark_accepted} != {values['modality_barcode_accepted_pairs']}")


def build_report_model(
    input_root: Path,
    library_name: str = "unknown library",
    pipeline_version: str = "unknown",
    run_metadata: Optional[dict[str, object]] = None,
) -> ReportModel:
    """Parse one completed output tree or a staged explicit-input directory."""
    files = find_files(input_root)
    by_name: dict[str, list[Path]] = defaultdict(list)
    for path in files:
        by_name[path.name].append(path)
    branches = discover_branches(files)
    warnings: list[str] = []

    marginal_stats: dict[tuple[str, str], tuple[dict, Path]] = {}
    cell_stats: dict[tuple[str, str], list[tuple[int, dict, Path]]] = defaultdict(list)
    retention_files: dict[str, list[dict[str, str]]] = {}
    gates: dict[tuple[str, str], tuple[list[dict[str, str]], Path]] = {}
    artifacts: dict[str, dict[str, object]] = {}
    star_logs: dict[str, tuple[dict, Path]] = {}
    summaries: dict[str, tuple[dict, Path]] = {}
    flagstats: dict[tuple[str, str, str], tuple[dict, Path]] = {}
    duplicates: dict[str, tuple[dict, Path]] = {}

    for path in files:
        name = path.name
        if match := re.match(r"^(.+)\.(rna|dna)_sample_barcode\.stats\.tsv$", name):
            marginal_stats[(match.group(1), match.group(2))] = (read_stats(path), path)
        elif match := re.match(r"^(.+)\.(rna|dna)_cell\.stats_L([123])\.tsv$", name):
            cell_stats[(match.group(1), match.group(2))].append((int(match.group(3)), read_stats(path), path))
        elif match := re.match(r"^(.+)\.dna_modality\.stats\.tsv$", name):
            marginal_stats[(match.group(1), "dna_modality")] = (read_stats(path), path)
        elif name.endswith((".rna_read_retention.tsv", ".dna_read_retention.tsv", ".rna_filter_retention.tsv", ".dna_alignment_retention.tsv")):
            retention_files[name] = read_tsv(path)
        elif match := re.match(r"^(.+)\.(rna|dna)_barcode_gates\.tsv$", name):
            gates[(match.group(1), match.group(2))] = (read_tsv(path), path)
        elif match := re.match(r"^(.+)\.dual_tag_artifact_filter\.summary\.tsv$", name):
            artifacts[match.group(1)] = read_artifact_summary(path)
        elif name.endswith(".Log.final.out"):
            star_logs[name.removesuffix(".Log.final.out")] = (read_star_log(path), path)
        elif name == "Summary.csv" and path.parent.name.endswith(".Solo.outGeneFull"):
            summaries[path.parent.name.removesuffix(".Solo.outGeneFull")] = (read_summary(path), path)
        elif match := re.match(r"^(rna|dna)\.(.+)\.(filtered_cells|aligned|markeddup|nodup)\.flagstat$", name):
            flagstats[(match.group(1), match.group(2), match.group(3))] = (read_flagstat(path), path)
        elif name.endswith(".DuplicateMetrics.txt"):
            duplicates[name.removesuffix(".DuplicateMetrics.txt")] = (read_duplicate_metrics(path), path)

    sidecar_composition = _parse_composition(files, branches, warnings)
    # The cumulative-gate sidecars remain authoritative for their own exact
    # denominator and are validated before sample composition is replaced by
    # the Tag counts/stats contract requested for the report visualization.
    _validate_composition_against_gates(sidecar_composition, gates)
    count_composition, count_contracts = _sample_composition_from_counts(files, branches)
    composition = [
        row
        for row in sidecar_composition
        if row["barcode_type"] != "sample_barcode"
        or (str(row["parent_sample"]), str(row["modality"]).lower()) not in count_contracts
    ]
    composition.extend(count_composition)
    legacy_sample_contracts = {
        (str(row["parent_sample"]), str(row["modality"]))
        for row in composition
        if row["barcode_type"] == "sample_barcode" and row.get("contract_level") == "legacy_group"
    }
    if legacy_sample_contracts:
        warnings.append(
            "Exact per-barcode composition unavailable for "
            + ", ".join(f"{parent} {modality}" for parent, modality in sorted(legacy_sample_contracts))
            + "; grouped legacy rows remain in the table and sample-composition plots are omitted."
        )
    retention: list[dict[str, object]] = []
    qc: list[dict[str, object]] = []
    complexity: list[dict[str, object]] = []
    common_done: set[tuple[str, str]] = set()
    artifact_done: set[str] = set()

    for branch in branches:
        raw, raw_source = _raw_pairs(marginal_stats, branch.parent, branch.modality)
        if raw is None:
            warnings.append(f"Raw input count unavailable for {branch.parent} {branch.modality}; retention path omitted.")
            continue
        common_key = (branch.parent, branch.modality)
        gate_contract = gates.get((branch.parent, branch.modality.lower()))
        gate_values = _gate_lookup(gate_contract[0], branch.modality.lower()) if gate_contract else None
        if common_key not in common_done:
            marginal = marginal_stats.get((branch.parent, branch.modality.lower()))
            if marginal:
                _collect_marginal_qc(qc, branch, marginal[0], marginal[1], "sample_barcode")
            modality_stats = marginal_stats.get((branch.parent, "dna_modality"))
            if modality_stats:
                _collect_marginal_qc(qc, branch, modality_stats[0], modality_stats[1], "modality_barcode")
            for layer, values, path in sorted(cell_stats.get((branch.parent, branch.modality.lower()), [])):
                matched = values.get(f"reads_with_L{layer}", {})
                _add_qc(qc, branch, f"cell_barcode_L{layer}_matched_pairs", matched.get("count"), "read_pairs", "raw input pairs; marginal layer match", path.name, "ALL")
                _add_qc(qc, branch, f"cell_barcode_L{layer}_match_pct", matched.get("percent"), "percent", "raw input pairs; marginal layer match", path.name, "ALL")
            if gate_values:
                for metric, value in gate_values.items():
                    _add_qc(qc, branch, metric, value, "read_pairs", "split-input pairs; exact cumulative gate", gate_contract[1].name, "ALL")
            common_done.add(common_key)

        stages: list[tuple] = []
        if branch.modality == "RNA":
            split_name = f"{branch.parent}.rna_read_retention.tsv"
            filter_name = f"{branch.split_id}.rna_filter_retention.tsv"
            split_rows = retention_files.get(split_name, [])
            filter_rows = retention_files.get(filter_name, [])
            split = _retention_lookup(split_rows)
            filt = _metric_lookup(filter_rows)
            if gate_values:
                split_input = split.get(("__all__", "__all__", "split_input_pairs"))
                if split_rows and split_input != gate_values["split_input_pairs"]:
                    raise ValueError(f"RNA split/gate input mismatch for {branch.parent}: {split_input} != {gate_values['split_input_pairs']}")
                stages.extend(
                    [
                        ("After paired trimming", gate_values["split_input_pairs"], "shared", gate_contract[1].name, "read pairs"),
                        ("L1–L2–L3 accepted", gate_values["ligation_barcode_accepted_pairs"], "shared", gate_contract[1].name, "read pairs"),
                        ("Sample barcode accepted — all required barcodes accepted", gate_values["sample_barcode_accepted_pairs"], "shared", gate_contract[1].name, "read pairs"),
                    ]
                )
                for label, metric in (("Properly paired", "paired_filter_pairs"), ("Canonical chromosomes", "canonical_pairs"), ("Called-cell final BAM", "called_cell_pairs")):
                    if metric in filt:
                        stages.append((label, filt[metric], "branch", by_name[filter_name][0].name, "primary R1 pair representatives"))
            else:
                warning = f"Exact RNA barcode gates unavailable for {branch.parent}; only provable coarse stages are shown."
                if warning not in warnings:
                    warnings.append(warning)
                split_input = split.get(("__all__", "__all__", "split_input_pairs"), split.get(("__all__", "__all__", "trimmed_input_pairs")))
                joint = split.get(("__all__", "__all__", "joint_barcode_accepted_pairs"))
                if split_input is not None:
                    stages.append(("After paired trimming", split_input, "shared", by_name[split_name][0].name, "read pairs"))
                if joint is not None:
                    stages.append(("Joint barcode accepted", joint, "shared", by_name[split_name][0].name, "read pairs"))
                star = star_logs.get(branch.split_id)
                if star:
                    star_input = integer(star[0].get("Number of input reads"))
                    mapped = (integer(star[0].get("Uniquely mapped reads number")) or 0) + (integer(star[0].get("Number of reads mapped to multiple loci")) or 0)
                    if star_input is not None:
                        stages.extend([("Entered STAR (coarse)", star_input, "branch", star[1].name, "read pairs"), ("STAR mapped primary (coarse)", mapped, "branch", star[1].name, "read pairs")])
                final = flagstats.get(("rna", branch.split_id, "filtered_cells"))
                if final and final[0].get("read1") is not None:
                    stages.append(("Called-cell final BAM", final[0]["read1"], "branch", final[1].name, "primary R1 pair representatives"))
            if stages:
                _add_path(retention, branch, raw, raw_source, stages)

            star = star_logs.get(branch.split_id)
            if star:
                for metric, source_key, unit, denominator in (
                    ("star_input_pairs", "Number of input reads", "read_pairs", "routed RNA branch"),
                    ("unique_genome_mapped_pairs", "Uniquely mapped reads number", "read_pairs", "STAR input pairs"),
                    ("unique_genome_mapping_pct", "Uniquely mapped reads %", "percent", "STAR input pairs"),
                ):
                    _add_qc(qc, branch, metric, star[0].get(source_key), unit, denominator, star[1].name)
            summary = summaries.get(branch.split_id)
            if summary:
                values, path = summary
                cells = integer(values.get("Estimated Number of Cells"))
                _add_qc(qc, branch, "estimated_cells", cells, "cells", "called STARsolo cells", path.name)
                saturation_fraction = numeric(values.get("Sequencing Saturation"))
                if saturation_fraction is not None and not 0 <= saturation_fraction <= 1:
                    raise ValueError(f"STARsolo sequencing saturation outside 0–1 in {path}")
                _add_qc(
                    qc,
                    branch,
                    "sequencing_saturation_pct",
                    saturation_fraction * 100 if saturation_fraction is not None else None,
                    "percent",
                    "STARsolo 2.7.11b: 1 - collapsed UMIs / reads assigned to a unique GeneFull feature",
                    path.name,
                )
                _add_qc(qc, branch, "median_umi_per_cell", values.get("Median UMI per Cell"), "UMIs_per_cell", "called STARsolo cells", path.name)
        else:
            split_name = f"{branch.parent}.dna_read_retention.tsv"
            align_name = f"{branch.split_id}.dna_alignment_retention.tsv"
            split_rows = retention_files.get(split_name, [])
            align_rows = retention_files.get(align_name, [])
            split = _retention_lookup(split_rows)
            align = _metric_lookup(align_rows)
            artifact = artifacts.get(branch.parent)
            if artifact and branch.parent not in artifact_done:
                values = [
                    ("dual_tag_filter_input_pairs", artifact["input_pairs"], "read_pairs", "pairs surviving paired trimming"),
                    ("dual_tag_filter_retained_pairs", artifact["retained_pairs"], "read_pairs", "dual-tag filter input pairs"),
                    ("dual_tag_filter_rejected_pairs", artifact["rejected_pairs"], "read_pairs", "dual-tag filter input pairs"),
                    ("dual_tag_filter_retained_pct", percentage(artifact["retained_pairs"], artifact["input_pairs"]), "percent", "dual-tag filter input pairs"),
                    ("dual_tag_filter_rejected_pct", percentage(artifact["rejected_pairs"], artifact["input_pairs"]), "percent", "dual-tag filter input pairs"),
                    ("dual_tag_filter_r1_signature_pairs", artifact["r1_with_signature"], "read_pairs", "dual-tag filter input pairs"),
                    ("dual_tag_filter_r2_signature_pairs", artifact["r2_with_signature"], "read_pairs", "dual-tag filter input pairs"),
                    ("dual_tag_filter_both_mates_signature_pairs", artifact["both_mates"], "read_pairs", "rejected dual-tag pairs"),
                    ("dual_tag_filter_r1_only_signature_pairs", artifact["r1_only"], "read_pairs", "rejected dual-tag pairs"),
                    ("dual_tag_filter_r2_only_signature_pairs", artifact["r2_only"], "read_pairs", "rejected dual-tag pairs"),
                    ("dual_tag_filter_cutadapt_version", artifact["cutadapt_version"], "version", "audit metadata"),
                    ("dual_tag_filter_signature_fasta_sha256", artifact["signature_fasta_sha256"], "sha256", "audit metadata"),
                ]
                for metric, value, unit, denominator in values:
                    _add_qc(qc, branch, metric, value, unit, denominator, artifact["path"], "ALL")
                artifact_done.add(branch.parent)
            if gate_values:
                split_input = split.get(("__all__", "__all__", "split_input_pairs"))
                if split_rows and split_input != gate_values["split_input_pairs"]:
                    raise ValueError(f"DNA split/gate input mismatch for {branch.parent}: {split_input} != {gate_values['split_input_pairs']}")
                if artifact:
                    if int(artifact["retained_pairs"]) != gate_values["split_input_pairs"]:
                        raise ValueError(f"Dual-tag retained/split input mismatch for {branch.parent}")
                    stages.extend([("After paired trimming", int(artifact["input_pairs"]), "shared", artifact["path"], "read pairs"), ("After DT artifact filter", int(artifact["retained_pairs"]), "shared", artifact["path"], "read pairs")])
                else:
                    stages.append(("After paired trimming", gate_values["split_input_pairs"], "shared", gate_contract[1].name, "read pairs"))
                stages.extend(
                    [
                        ("L1–L2–L3 accepted", gate_values["ligation_barcode_accepted_pairs"], "shared", gate_contract[1].name, "read pairs"),
                        ("Sample barcode accepted", gate_values["sample_barcode_accepted_pairs"], "shared", gate_contract[1].name, "read pairs"),
                        ("MO barcode accepted — all required barcodes accepted", gate_values["modality_barcode_accepted_pairs"], "shared", gate_contract[1].name, "read pairs"),
                    ]
                )
                routed = split.get((branch.group, branch.branch, "routed_branch_pairs"))
                if routed is not None:
                    stages.append(("Mark-specific routing branch", routed, "branch", by_name[split_name][0].name, "read pairs"))
                if "proper_pair_primary_pairs" in align:
                    stages.append(("Properly paired after blacklist", align["proper_pair_primary_pairs"], "branch", by_name[align_name][0].name, "primary R1 pair representatives"))
            else:
                warning = f"Exact DNA barcode gates unavailable for {branch.parent}; only provable coarse stages are shown."
                if warning not in warnings:
                    warnings.append(warning)
                split_input = split.get(("__all__", "__all__", "split_input_pairs"), split.get(("__all__", "__all__", "trimmed_input_pairs")))
                joint = split.get(("__all__", "__all__", "joint_barcode_accepted_pairs"))
                if artifact and split_input is not None:
                    if int(artifact["retained_pairs"]) != split_input:
                        raise ValueError(f"Legacy dual-tag retained/split input mismatch for {branch.parent}")
                    stages.extend([("After paired trimming", int(artifact["input_pairs"]), "shared", artifact["path"], "read pairs"), ("After DT artifact filter", int(artifact["retained_pairs"]), "shared", artifact["path"], "read pairs")])
                elif split_input is not None:
                    stages.append(("After paired trimming", split_input, "shared", by_name[split_name][0].name, "read pairs"))
                if joint is not None:
                    stages.append(("Joint barcode accepted", joint, "shared", by_name[split_name][0].name, "read pairs"))
                duplicate_entry = duplicates.get(branch.split_id)
                if duplicate_entry and duplicate_entry[0].get("read_pairs_examined") is not None:
                    stages.append(("Proper-pair + blacklist (coarse)", duplicate_entry[0]["read_pairs_examined"], "branch", duplicate_entry[1].name, "read pairs"))
            canonical = flagstats.get(("dna", branch.split_id, "markeddup")) or flagstats.get(("dna", branch.split_id, "aligned"))
            nodup = flagstats.get(("dna", branch.split_id, "nodup"))
            if canonical and canonical[0].get("read1") is not None:
                stages.append(("Canonical marked-duplicate BAM", canonical[0]["read1"], "branch", canonical[1].name, "primary R1 pair representatives"))
            if nodup and nodup[0].get("read1") is not None:
                stages.append(("Canonical NoDup final", nodup[0]["read1"], "branch", nodup[1].name, "primary R1 pair representatives"))
            if stages:
                _add_path(retention, branch, raw, raw_source, stages)

            duplicate_entry = duplicates.get(branch.split_id)
            if duplicate_entry and "read_pairs_examined" in duplicate_entry[0]:
                values, path = duplicate_entry
                examined = integer(values.get("read_pairs_examined"))
                duplicate_pairs = integer(values.get("read_pair_duplicates"))
                optical = integer(values.get("read_pair_optical_duplicates"))
                if None in {examined, duplicate_pairs, optical}:
                    raise ValueError(f"Incomplete Picard duplication metrics in {path}")
                unique = examined - duplicate_pairs
                pcr_library = duplicate_pairs - optical
                if min(unique, pcr_library, optical) < 0 or unique + pcr_library + optical != examined:
                    raise ValueError(f"Picard duplicate components do not reconcile in {path}")
                summary_values = [
                    ("read_pairs_examined", examined, "read_pairs"),
                    ("read_pair_duplicates", duplicate_pairs, "read_pairs"),
                    ("read_pair_optical_duplicates", optical, "read_pairs"),
                    ("observed_unique_pairs", unique, "read_pairs"),
                    ("pcr_or_library_duplicate_pairs", pcr_library, "read_pairs"),
                    ("percent_duplication", values.get("percent_duplication"), "percent"),
                    ("estimated_library_size", values.get("estimated_library_size"), "fragments"),
                ]
                for metric, value, unit in summary_values:
                    _add_qc(qc, branch, metric, value, unit, "Picard examined pairs" if metric != "estimated_library_size" else "Picard estimate", path.name)
                points = values.get("roi") or [None]
                for point in points:
                    multiplier = point["yield_multiplier"] if point else None
                    complexity.append(
                        {
                            "run": branch.run, "parent_sample": branch.parent, "branch": branch.branch,
                            "read_pairs_examined": examined, "observed_unique_pairs": unique,
                            "pcr_or_library_duplicate_pairs": pcr_library, "optical_duplicate_pairs": optical,
                            "observed_unique_fraction": fraction(unique, examined),
                            "pcr_or_library_duplicate_fraction": fraction(pcr_library, examined),
                            "optical_duplicate_fraction": fraction(optical, examined),
                            "percent_duplication": values.get("percent_duplication"),
                            "estimated_library_size": values.get("estimated_library_size"),
                            "roi_depth_multiplier": point["depth_multiplier"] if point else None,
                            "roi_yield_multiplier": multiplier,
                            "roi_estimated_unique_pairs": unique * multiplier if multiplier is not None else None,
                            "source": path.name,
                        }
                    )

    if not composition:
        warnings.append("Exact barcode-composition counters are unavailable; composition sections are omitted.")
    raw_duplicates = [values for values, _ in duplicates.values()]
    raw_artifacts = [{key: value for key, value in artifact.items() if key != "path"} for artifact in artifacts.values()]
    return ReportModel(
        library_name=library_name or "unknown library",
        pipeline_version=pipeline_version or "unknown",
        generated_at=datetime.now(timezone.utc).isoformat(),
        branches=branches,
        retention=retention,
        qc_metrics=qc,
        barcode_composition=composition,
        library_complexity=complexity,
        insights=[],
        warnings=list(dict.fromkeys(warnings)),
        run_metadata=run_metadata or {},
        raw={"duplicate_metrics": raw_duplicates, "dual_tag_artifact_filters": raw_artifacts},
    )


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            normalized = dict(row)
            for key in ("retained_prev_pct", "cumulative_raw_pct", "percentage", "observed_unique_fraction", "pcr_or_library_duplicate_fraction", "optical_duplicate_fraction", "percent_duplication"):
                if key in normalized and normalized[key] is not None:
                    normalized[key] = f"{float(normalized[key]):.6f}"
            writer.writerow(normalized)


def write_model_tables(model: ReportModel, output_dir: Path) -> dict[str, Path]:
    outputs = {
        "read_retention": output_dir / "read_retention.tsv",
        "qc_metrics": output_dir / "qc_metrics.tsv",
        "barcode_composition": output_dir / "barcode_composition.tsv",
        "library_complexity": output_dir / "library_complexity.tsv",
    }
    write_tsv(outputs["read_retention"], RETENTION_FIELDS, model.retention)
    write_tsv(outputs["qc_metrics"], QC_FIELDS, model.qc_metrics)
    write_tsv(outputs["barcode_composition"], COMPOSITION_FIELDS, model.barcode_composition)
    write_tsv(outputs["library_complexity"], COMPLEXITY_FIELDS, model.library_complexity)
    return outputs


def write_model_json(model: ReportModel, path: Path) -> None:
    path.write_text(json.dumps(model.to_dict(), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
