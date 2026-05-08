#!/usr/bin/env python3
"""Build RNA/DNA sequencing-efficiency tables and plots from TrESFlow outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):[A-Za-z]:(.*)$")

RNA_TAG_STAGES = [
    ("valid_sb", "valid sample barcode"),
    ("valid_l1", "valid L1 barcode"),
    ("valid_l2", "valid L2 barcode"),
    ("valid_l3", "valid L3 barcode"),
    ("valid_cb", "valid full cell barcode"),
    ("umi_present", "UMI present"),
]

DNA_TAG_STAGES = [
    ("valid_sb", "valid sample barcode"),
    ("valid_l1", "valid L1 barcode"),
    ("valid_l2", "valid L2 barcode"),
    ("valid_l3", "valid L3 barcode"),
    ("valid_cb", "valid full cell barcode"),
    ("valid_mo", "valid modality barcode"),
]

RNA_BAM_STAGES = [
    ("aligned", "aligned reads from filtered_cells BAM"),
    ("gx_assigned", "gene-assigned reads"),
]

DNA_BAM_STAGES = [
    ("pre_dedup_aligned", "aligned reads before duplicate removal"),
    ("post_dedup_retained", "aligned reads after duplicate removal"),
]

UPSET_LABELS = {
    "valid_sb": "valid SB",
    "valid_l1": "valid L1",
    "valid_l2": "valid L2",
    "valid_l3": "valid L3",
    "valid_cb": "valid CB",
    "umi_present": "UMI present",
    "valid_mo": "valid MO",
    "aligned": "aligned",
    "gx_assigned": "GX assigned",
    "pre_dedup_aligned": "pre-dedup aligned",
    "post_dedup_retained": "post-dedup retained",
}


@dataclass
class WarningRecord:
    modality: str
    unit: str
    message: str


@dataclass
class UnitReport:
    name: str
    modality: str
    level: str
    tagged_ids: Set[str] = field(default_factory=set)
    stage_sets: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    stage_available: Set[str] = field(default_factory=set)
    warnings: List[str] = field(default_factory=list)

    def add_warning(self, message: str) -> None:
        if message not in self.warnings:
            self.warnings.append(message)

    def add_tag_record(self, read_id: str, tags: Dict[str, str]) -> None:
        self.tagged_ids.add(read_id)
        checks = {
            "valid_sb": is_present(tags.get("SB")),
            "valid_l1": is_present(tags.get("L1")),
            "valid_l2": is_present(tags.get("L2")),
            "valid_l3": is_present(tags.get("L3")),
            "valid_cb": is_present(tags.get("CB")),
            "umi_present": is_present(tags.get("UM")),
            "valid_mo": is_present(tags.get("MO")),
        }

        for stage, passed in checks.items():
            if passed:
                self.stage_sets[stage].add(read_id)

    def add_stage_ids(self, stage: str, read_ids: Set[str]) -> None:
        self.stage_available.add(stage)
        self.stage_sets[stage].update(read_ids)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--rna-tag-records", nargs="*", default=[], type=Path)
    parser.add_argument("--rna-filtered-bams", nargs="*", default=[], type=Path)
    parser.add_argument("--dna-tag-records", nargs="*", default=[], type=Path)
    parser.add_argument("--dna-markeddup-bams", nargs="*", default=[], type=Path)
    parser.add_argument("--dna-nodup-bams", nargs="*", default=[], type=Path)
    parser.add_argument("--sb-group-maps", nargs="*", default=[], type=Path)
    parser.add_argument("--dna-mo-maps", nargs="*", default=[], type=Path)
    return parser.parse_args()


def is_present(value: Optional[str]) -> bool:
    if value is None:
        return False
    cleaned = str(value).strip()
    return cleaned not in {"", "-", "NoMatch", "nomatch", "NA", "N/A", "None", "null"}


def read_pairs(read_records: Optional[int]) -> Optional[float]:
    if read_records is None:
        return None
    return read_records / 2.0


def format_number(value: Optional[float]) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_tag_record_line(line: str) -> Tuple[Optional[str], Dict[str, str]]:
    fields = line.rstrip("\n").split("\t")
    if not fields or not fields[0].strip():
        return None, {}
    read_id = fields[0].strip()
    tags: Dict[str, str] = {}
    for field in fields[1:]:
        for token in field.replace("\t", " ").split():
            match = TAG_RE.match(token)
            if match:
                tags[match.group(1)] = match.group(2)
    return read_id, tags


def sample_from_tag_record(path: Path, modality: str) -> str:
    name = path.name
    suffixes = [
        f".{modality}_tag_records.tsv.gz",
        f".{modality}_tag_records.tsv",
        ".tag_records.tsv.gz",
        ".tag_records.tsv",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name.split(".")[0]


def load_sb_group_maps(paths: Sequence[Path]) -> Tuple[Dict[str, Dict[str, str]], Dict[str, List[str]]]:
    sb_to_group: Dict[str, Dict[str, str]] = defaultdict(dict)
    groups: Dict[str, List[str]] = defaultdict(list)
    seen_groups: Dict[str, Set[str]] = defaultdict(set)

    for path in paths:
        if not path.exists():
            continue
        with open(path, "rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample = (row.get("sample") or "").strip()
                group = (row.get("sb_group") or "").strip()
                sb_bc = (row.get("sb_bc") or "").strip()
                if not sample or not group or not sb_bc:
                    continue
                sb_to_group[sample][sb_bc] = group
                if group not in seen_groups[sample]:
                    groups[sample].append(group)
                    seen_groups[sample].add(group)

    return dict(sb_to_group), dict(groups)


def split_sb_group_map_paths(paths: Sequence[Path]) -> Tuple[List[Path], List[Path]]:
    rna_paths = [path for path in paths if path.name.startswith("rna_")]
    dna_paths = [path for path in paths if path.name.startswith("dna_")]
    generic_paths = [path for path in paths if path not in rna_paths and path not in dna_paths]
    return (rna_paths or generic_paths, dna_paths or generic_paths)


def load_dna_mo_maps(paths: Sequence[Path]) -> Dict[str, Dict[Tuple[str, str], str]]:
    mapping: Dict[str, Dict[Tuple[str, str], str]] = defaultdict(dict)

    for path in paths:
        if not path.exists():
            continue
        with open(path, "rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample = (row.get("sample") or "").strip()
                group = (row.get("sb_group") or "").strip()
                mark = (row.get("mark") or "").strip()
                mo_bc = (row.get("mo_bc") or "").strip()
                if not sample or not group or not mark or not mo_bc:
                    continue
                mapping[sample][(group, mo_bc)] = mark

    return dict(mapping)


def normalize_sb_drop_first(value: str) -> str:
    return value[1:] if len(value) > 1 else value


def resolve_group(sample: str, sb_value: Optional[str], sb_maps: Dict[str, Dict[str, str]]) -> Optional[str]:
    if not is_present(sb_value):
        return None
    sample_map = sb_maps.get(sample, {})
    if sb_value in sample_map:
        return sample_map[sb_value]
    trimmed = normalize_sb_drop_first(str(sb_value))
    return sample_map.get(trimmed)


def resolve_mark(
    sample: str,
    group: Optional[str],
    mo_value: Optional[str],
    mo_maps: Dict[str, Dict[Tuple[str, str], str]],
) -> Optional[str]:
    if not group or not is_present(mo_value):
        return None
    sample_map = mo_maps.get(sample, {})
    return sample_map.get((group, str(mo_value)))


def unit_key(modality: str, name: str) -> Tuple[str, str]:
    return modality, name


def get_unit(units: Dict[Tuple[str, str], UnitReport], modality: str, name: str, level: str) -> UnitReport:
    key = unit_key(modality, name)
    if key not in units:
        units[key] = UnitReport(name=name, modality=modality, level=level)
    return units[key]


def add_tag_records(
    units: Dict[Tuple[str, str], UnitReport],
    modality: str,
    path: Path,
    sb_maps: Dict[str, Dict[str, str]],
    mo_maps: Dict[str, Dict[Tuple[str, str], str]],
    warnings: List[WarningRecord],
) -> None:
    sample = sample_from_tag_record(path, modality)
    sample_unit = get_unit(units, modality, sample, "sample")

    try:
        with open_maybe_gzip(path) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                read_id, tags = parse_tag_record_line(raw_line)
                if read_id is None:
                    continue

                sample_unit.add_tag_record(read_id, tags)
                group = resolve_group(sample, tags.get("SB"), sb_maps)
                if group:
                    group_unit = get_unit(units, modality, f"{sample}_{group}", "sample_group")
                    group_unit.add_tag_record(read_id, tags)
                    if modality == "dna":
                        mark = resolve_mark(sample, group, tags.get("MO"), mo_maps)
                        if mark:
                            mark_unit = get_unit(units, modality, f"{sample}_{group}_{mark}", "sample_group_mark")
                            mark_unit.add_tag_record(read_id, tags)
    except Exception as exc:  # pragma: no cover - exercised by integration tests
        message = f"Could not parse tag records {path}: {exc}"
        sample_unit.add_warning(message)
        warnings.append(WarningRecord(modality, sample, message))


def parse_split_name(split_name: str, sample_ids: Iterable[str], modality: str) -> Tuple[str, Optional[str], Optional[str]]:
    for sample in sorted(set(sample_ids), key=len, reverse=True):
        if split_name == sample:
            return sample, None, None
        prefix = f"{sample}_"
        if split_name.startswith(prefix):
            suffix = split_name[len(prefix) :]
            if modality == "rna":
                return sample, suffix or None, None
            parts = suffix.split("_")
            if len(parts) >= 2:
                return sample, parts[0], "_".join(parts[1:])
            return sample, parts[0] if parts else None, None

    tokens = split_name.split("_")
    if modality == "dna" and len(tokens) >= 3:
        return "_".join(tokens[:-2]), tokens[-2], tokens[-1]
    if modality == "rna" and len(tokens) >= 2:
        return "_".join(tokens[:-1]), tokens[-1], None
    return split_name, None, None


def split_from_bam(path: Path, suffix: str) -> str:
    name = path.name
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return name


def safe_bam_ids(path: Path, gx_only: bool = False) -> Tuple[Optional[Set[str]], Optional[str]]:
    try:
        import pysam  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on runtime env
        return None, f"pysam is unavailable; skipping BAM-derived stage for {path}: {exc}"

    read_ids: Set[str] = set()
    try:
        with pysam.AlignmentFile(str(path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                if gx_only:
                    if not read.has_tag("GX"):
                        continue
                    gx_value = read.get_tag("GX")
                    if not is_present(str(gx_value)):
                        continue
                read_ids.add(read.query_name)
    except Exception as exc:
        return None, f"Could not read BAM-derived stage from {path}; skipping stage: {exc}"

    return read_ids, None


def add_bam_stage(
    units: Dict[Tuple[str, str], UnitReport],
    modality: str,
    path: Path,
    suffix: str,
    stage: str,
    sample_ids: Iterable[str],
    warnings: List[WarningRecord],
    gx_only: bool = False,
) -> None:
    split_name = split_from_bam(path, suffix)
    sample, group, mark = parse_split_name(split_name, sample_ids, modality)
    target_units = [get_unit(units, modality, sample, "sample")]
    if group:
        target_units.append(get_unit(units, modality, f"{sample}_{group}", "sample_group"))
    if modality == "dna" and group and mark:
        target_units.append(get_unit(units, modality, f"{sample}_{group}_{mark}", "sample_group_mark"))

    read_ids, warning = safe_bam_ids(path, gx_only=gx_only)
    if warning:
        for unit in target_units:
            unit.add_warning(warning)
            warnings.append(WarningRecord(modality, unit.name, warning))
        return

    assert read_ids is not None
    for unit in target_units:
        unit.add_stage_ids(stage, read_ids)


def stage_defs_for_modality(modality: str) -> List[Tuple[str, str]]:
    if modality == "rna":
        return [("total_tagged", "total tagged records")] + RNA_TAG_STAGES + RNA_BAM_STAGES + [
            ("final_passing", "final passing reads")
        ]
    return [("total_tagged", "total tagged records")] + DNA_TAG_STAGES + DNA_BAM_STAGES + [
        ("final_passing", "final passing reads")
    ]


def final_stage_name(modality: str) -> str:
    return "gx_assigned" if modality == "rna" else "post_dedup_retained"


def available_tag_stage(unit: UnitReport, stage: str) -> bool:
    return bool(unit.tagged_ids) and stage in unit.stage_sets


def make_stage_rows(unit: UnitReport) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    active_ids: Optional[Set[str]] = set(unit.tagged_ids) if unit.tagged_ids else None

    for order, (stage, label) in enumerate(stage_defs_for_modality(unit.modality), start=1):
        count: Optional[int]
        available = True
        note = ""

        if stage == "total_tagged":
            count = len(unit.tagged_ids) if unit.tagged_ids else None
            available = bool(unit.tagged_ids)
        elif stage == "final_passing":
            final_stage = final_stage_name(unit.modality)
            if final_stage in unit.stage_sets and final_stage in unit.stage_available:
                count = len(unit.stage_sets[final_stage])
            elif final_stage in unit.stage_sets and active_ids is not None:
                count = len(active_ids & unit.stage_sets[final_stage])
            else:
                count = None
                available = False
                note = "Final BAM-derived stage unavailable"
        else:
            is_bam_stage = stage in {name for name, _ in RNA_BAM_STAGES + DNA_BAM_STAGES}
            if is_bam_stage and stage not in unit.stage_available:
                count = None
                available = False
                note = "Optional BAM-derived stage unavailable"
            elif stage not in unit.stage_sets:
                count = None
                available = False
                note = "No reads passed this stage"
            elif active_ids is not None:
                stage_ids = unit.stage_sets[stage]
                intersected = active_ids & stage_ids
                if is_bam_stage and active_ids and stage_ids and not intersected:
                    count = len(stage_ids)
                    active_ids = set(stage_ids)
                    note = "Read-name intersection unavailable; using BAM stage count"
                else:
                    count = len(intersected)
                    active_ids = intersected
            else:
                count = len(unit.stage_sets[stage])
                active_ids = set(unit.stage_sets[stage])

        if not available and count is None:
            active_ids = None if stage in {name for name, _ in RNA_BAM_STAGES + DNA_BAM_STAGES} else active_ids

        rows.append(
            {
                "unit": unit.name,
                "unit_level": unit.level,
                "modality": unit.modality,
                "stage_order": order,
                "stage": stage,
                "stage_label": label,
                "read_records": count,
                "read_pairs": read_pairs(count),
                "available": str(bool(available)).lower(),
                "note": note,
            }
        )

    return rows


def write_table(path: Path, rows: List[Dict[str, object]]) -> None:
    fieldnames = [
        "unit",
        "unit_level",
        "modality",
        "stage_order",
        "stage",
        "stage_label",
        "read_records",
        "read_pairs",
        "available",
        "note",
    ]
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            formatted["read_records"] = format_number(formatted["read_records"])  # type: ignore[arg-type]
            formatted["read_pairs"] = format_number(formatted["read_pairs"])  # type: ignore[arg-type]
            writer.writerow(formatted)


def available_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    return [row for row in rows if row["available"] == "true" and row["read_records"] is not None]


def write_placeholder_html(path: Path, title: str, message: str) -> None:
    path.write_text(
        "<html><head><meta charset=\"utf-8\"><title>"
        + html.escape(title)
        + "</title></head><body><h1>"
        + html.escape(title)
        + "</h1><p>"
        + html.escape(message)
        + "</p></body></html>\n",
        encoding="utf-8",
    )


def write_placeholder_pdf(path: Path, title: str, message: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis("off")
        ax.text(0.02, 0.8, title, fontsize=14, weight="bold", transform=ax.transAxes)
        ax.text(0.02, 0.55, message, fontsize=10, transform=ax.transAxes, wrap=True)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
    except Exception:
        path.write_text(f"{title}\n\n{message}\n", encoding="utf-8")


def write_sankey(unit: UnitReport, rows: List[Dict[str, object]], html_path: Path, pdf_path: Path) -> None:
    data_rows = available_rows(rows)
    title = f"{unit.name} {unit.modality.upper()} sequencing efficiency"
    if len(data_rows) < 2:
        message = "Not enough available stages to draw a sequencing-efficiency Sankey plot."
        write_placeholder_html(html_path, title, message)
        write_placeholder_pdf(pdf_path, title, message)
        return

    labels = [str(row["stage_label"]) for row in data_rows]
    counts = [int(row["read_records"]) for row in data_rows]  # type: ignore[arg-type]

    try:
        import plotly.graph_objects as go
        from plotly.offline import plot

        node_labels = list(labels)
        sources: List[int] = []
        targets: List[int] = []
        values: List[int] = []

        for idx in range(len(counts) - 1):
            current = counts[idx]
            nxt = counts[idx + 1]
            retained = max(min(current, nxt), 0)
            lost = max(current - retained, 0)
            sources.append(idx)
            targets.append(idx + 1)
            values.append(retained)
            if lost > 0:
                loss_idx = len(node_labels)
                node_labels.append(f"lost after {labels[idx]}")
                sources.append(idx)
                targets.append(loss_idx)
                values.append(lost)

        figure = go.Figure(
            data=[
                go.Sankey(
                    arrangement="snap",
                    node={"label": node_labels, "pad": 18, "thickness": 14},
                    link={"source": sources, "target": targets, "value": values},
                )
            ]
        )
        figure.update_layout(title_text=title, font_size=11)
        plot(figure, filename=str(html_path), auto_open=False, include_plotlyjs="cdn")
        try:
            figure.write_image(str(pdf_path))
            return
        except Exception:
            pass
    except Exception:
        write_placeholder_html(html_path, title, "Plotly is unavailable; see the PDF funnel plot.")

    write_funnel_pdf(pdf_path, title, labels, counts)


def write_funnel_pdf(path: Path, title: str, labels: List[str], counts: List[int]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        height = max(3.0, 0.42 * len(labels) + 1.5)
        fig, ax = plt.subplots(figsize=(9, height))
        y_labels = labels[::-1]
        y_counts = counts[::-1]
        ax.barh(y_labels, y_counts, color="#4C78A8")
        ax.set_xlabel("read records")
        ax.set_title(title)
        for idx, value in enumerate(y_counts):
            ax.text(value, idx, f" {value}", va="center", fontsize=8)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
    except Exception as exc:
        path.write_text(f"{title}\n\nUnable to render PDF plot: {exc}\n", encoding="utf-8")


def write_upset(unit: UnitReport, html_path: Path, pdf_path: Path) -> None:
    title = f"{unit.name} {unit.modality.upper()} sequencing-efficiency intersections"
    stage_order = [stage for stage, _ in (RNA_TAG_STAGES + RNA_BAM_STAGES if unit.modality == "rna" else DNA_TAG_STAGES + DNA_BAM_STAGES)]
    usable_stages = [stage for stage in stage_order if unit.stage_sets.get(stage)]

    if len(usable_stages) < 2:
        message = "Not enough available stage sets to draw an UpSet plot."
        write_placeholder_html(html_path, title, message)
        write_placeholder_pdf(pdf_path, title, message)
        return

    try:
        import pandas as pd

        universe = sorted(set().union(*(unit.stage_sets[stage] for stage in usable_stages)))
        if not universe:
            raise ValueError("No read identifiers available for UpSet intersections")

        labels = [UPSET_LABELS.get(stage, stage) for stage in usable_stages]
        data = {
            label: [read_id in unit.stage_sets[stage] for read_id in universe]
            for label, stage in zip(labels, usable_stages)
        }
        frame = pd.DataFrame(data)
        summary = frame.groupby(labels).size().rename("read_records").reset_index()
        summary["read_pairs"] = summary["read_records"] / 2.0
        summary.to_html(html_path, index=False)

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from upsetplot import UpSet

        upset_data = frame.assign(read_records=1).groupby(labels)["read_records"].sum()
        fig = plt.figure(figsize=(max(8, len(labels) * 1.2), 5))
        UpSet(upset_data, subset_size="count", show_counts=True).plot(fig=fig)
        fig.suptitle(title)
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
    except Exception as exc:
        message = f"UpSet plot skipped because exact intersections could not be rendered: {exc}"
        unit.add_warning(message)
        write_placeholder_html(html_path, title, message)
        write_placeholder_pdf(pdf_path, title, message)


def write_warnings(path: Path, warnings: Sequence[WarningRecord]) -> None:
    seen: Set[Tuple[str, str, str]] = set()
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["level", "modality", "unit", "message"])
        for warning in warnings:
            key = (warning.modality, warning.unit, warning.message)
            if key in seen:
                continue
            seen.add(key)
            writer.writerow(["WARNING", warning.modality, warning.unit, warning.message])
            print(f"WARNING [{warning.modality}:{warning.unit}] {warning.message}", file=sys.stderr)


def write_combined_summary(outdir: Path, units: Dict[Tuple[str, str], UnitReport]) -> None:
    samples = sorted({unit.name for unit in units.values() if unit.level == "sample"})
    for sample in samples:
        rows = []
        for modality in ("rna", "dna"):
            unit = units.get(unit_key(modality, sample))
            if unit is None:
                continue
            stage_rows = make_stage_rows(unit)
            total = next((row for row in stage_rows if row["stage"] == "total_tagged"), None)
            final = next((row for row in stage_rows if row["stage"] == "final_passing"), None)
            total_count = total["read_records"] if total else None
            final_count = final["read_records"] if final else None
            retention = ""
            if isinstance(total_count, int) and total_count > 0 and isinstance(final_count, int):
                retention = format_number(final_count / total_count)
            rows.append(
                {
                    "sample": sample,
                    "modality": modality,
                    "total_tagged_read_records": total_count,
                    "total_tagged_read_pairs": read_pairs(total_count if isinstance(total_count, int) else None),
                    "final_passing_read_records": final_count,
                    "final_passing_read_pairs": read_pairs(final_count if isinstance(final_count, int) else None),
                    "final_retention_fraction": retention,
                }
            )

        if not rows:
            continue

        table_path = outdir / f"{sample}.combined_sequencing_efficiency.tsv"
        with open(table_path, "wt", encoding="utf-8", newline="") as handle:
            fieldnames = [
                "sample",
                "modality",
                "total_tagged_read_records",
                "total_tagged_read_pairs",
                "final_passing_read_records",
                "final_passing_read_pairs",
                "final_retention_fraction",
            ]
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                formatted = dict(row)
                for key in [
                    "total_tagged_read_records",
                    "total_tagged_read_pairs",
                    "final_passing_read_records",
                    "final_passing_read_pairs",
                ]:
                    formatted[key] = format_number(formatted[key])  # type: ignore[arg-type]
                writer.writerow(formatted)

        write_combined_html(outdir / f"{sample}.combined_sequencing_efficiency.html", sample, rows)
        write_combined_pdf(outdir / f"{sample}.combined_sequencing_efficiency.pdf", sample, rows)


def write_combined_html(path: Path, sample: str, rows: List[Dict[str, object]]) -> None:
    header = "".join(f"<th>{html.escape(key)}</th>" for key in rows[0].keys())
    body = []
    for row in rows:
        body.append("".join(f"<td>{html.escape(format_number(value) if isinstance(value, float) else str(value or ''))}</td>" for value in row.values()))
    path.write_text(
        "<html><head><meta charset=\"utf-8\"><title>"
        + html.escape(sample)
        + " combined sequencing efficiency</title></head><body><table><thead><tr>"
        + header
        + "</tr></thead><tbody>"
        + "".join(f"<tr>{cells}</tr>" for cells in body)
        + "</tbody></table></body></html>\n",
        encoding="utf-8",
    )


def write_combined_pdf(path: Path, sample: str, rows: List[Dict[str, object]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        labels = [str(row["modality"]).upper() for row in rows]
        totals = [int(row["total_tagged_read_records"] or 0) for row in rows]
        finals = [int(row["final_passing_read_records"] or 0) for row in rows]

        fig, ax = plt.subplots(figsize=(6, 4))
        x_values = range(len(labels))
        ax.bar([x - 0.18 for x in x_values], totals, width=0.36, label="total tagged")
        ax.bar([x + 0.18 for x in x_values], finals, width=0.36, label="final passing")
        ax.set_xticks(list(x_values), labels)
        ax.set_ylabel("read records")
        ax.set_title(f"{sample} combined sequencing efficiency")
        ax.legend()
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
    except Exception as exc:
        path.write_text(f"{sample} combined sequencing efficiency\n\nUnable to render PDF: {exc}\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rna_sb_paths, dna_sb_paths = split_sb_group_map_paths(args.sb_group_maps)
    rna_sb_maps, rna_groups_by_sample = load_sb_group_maps(rna_sb_paths)
    dna_sb_maps, dna_groups_by_sample = load_sb_group_maps(dna_sb_paths)
    mo_maps = load_dna_mo_maps(args.dna_mo_maps)
    units: Dict[Tuple[str, str], UnitReport] = {}
    warnings: List[WarningRecord] = []

    for sample, groups in rna_groups_by_sample.items():
        for group in groups:
            get_unit(units, "rna", f"{sample}_{group}", "sample_group")
    for sample, groups in dna_groups_by_sample.items():
        for group in groups:
            get_unit(units, "dna", f"{sample}_{group}", "sample_group")

    for path in args.rna_tag_records:
        add_tag_records(units, "rna", path, rna_sb_maps, mo_maps, warnings)
    for path in args.dna_tag_records:
        add_tag_records(units, "dna", path, dna_sb_maps, mo_maps, warnings)

    sample_ids = {sample_from_tag_record(path, "rna") for path in args.rna_tag_records}
    sample_ids.update(sample_from_tag_record(path, "dna") for path in args.dna_tag_records)
    sample_ids.update(rna_sb_maps.keys())
    sample_ids.update(dna_sb_maps.keys())

    for path in args.rna_filtered_bams:
        add_bam_stage(units, "rna", path, ".filtered_cells.bam", "aligned", sample_ids, warnings)
        add_bam_stage(units, "rna", path, ".filtered_cells.bam", "gx_assigned", sample_ids, warnings, gx_only=True)
    for path in args.dna_markeddup_bams:
        add_bam_stage(units, "dna", path, "_MarkedDup.bam", "pre_dedup_aligned", sample_ids, warnings)
    for path in args.dna_nodup_bams:
        add_bam_stage(units, "dna", path, "_NoDup.bam", "post_dedup_retained", sample_ids, warnings)

    all_unit_warnings: List[WarningRecord] = list(warnings)
    for unit in sorted(units.values(), key=lambda item: (item.modality, item.name)):
        if not unit.tagged_ids and not any(unit.stage_sets.values()):
            continue
        rows = make_stage_rows(unit)
        for row in rows:
            note = str(row.get("note") or "")
            if row.get("available") == "false" and (
                note.startswith("Optional BAM-derived stage unavailable")
                or note.startswith("Final BAM-derived stage unavailable")
            ):
                unit.add_warning(f"{row['stage_label']}: {note}")
        prefix = args.outdir / f"{unit.name}.{unit.modality}_sequencing_efficiency"
        write_table(Path(f"{prefix}.tsv"), rows)
        write_sankey(unit, rows, Path(f"{prefix}.sankey.html"), Path(f"{prefix}.sankey.pdf"))
        write_upset(unit, Path(f"{prefix}.upset.html"), Path(f"{prefix}.upset.pdf"))
        for message in unit.warnings:
            all_unit_warnings.append(WarningRecord(unit.modality, unit.name, message))

    write_combined_summary(args.outdir, units)
    write_warnings(args.outdir / "sequencing_efficiency.warnings.tsv", all_unit_warnings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
