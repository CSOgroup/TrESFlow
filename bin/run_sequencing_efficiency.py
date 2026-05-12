#!/usr/bin/env python3
"""Build RNA/DNA sequencing-efficiency UpSet PDFs from TrESFlow outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):[A-Za-z]:(.*)$")

COMMON_CATEGORIES = [
    ("reads", "Reads +"),
    ("sample", "Sample +"),
    ("ligation", "Ligation +"),
    ("cb", "CB +"),
    ("cb100", "CB>100 +"),
]

RNA_CATEGORIES = COMMON_CATEGORIES + [
    ("umi", "UMI +"),
    ("mapped", "Mapped +"),
    ("gx", "GX +"),
]

DNA_CATEGORIES = COMMON_CATEGORIES + [
    ("modality", "Modality +"),
    ("mapped", "Mapped +"),
    ("unique", "Unique +"),
]


@dataclass(frozen=True)
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
    warnings: Set[str] = field(default_factory=set)

    def add_warning(self, message: str) -> None:
        self.warnings.add(message)

    def register_tag_categories(self) -> None:
        categories = RNA_CATEGORIES if self.modality == "rna" else DNA_CATEGORIES
        for stage, _label in categories:
            if stage in {"mapped", "gx", "unique", "cb100"}:
                continue
            self.stage_available.add(stage)

    def add_tag_record(self, read_id: str, tags: Dict[str, str]) -> None:
        self.register_tag_categories()
        self.tagged_ids.add(read_id)
        self.stage_sets["reads"].add(read_id)

        if is_present(tags.get("SB")):
            self.stage_sets["sample"].add(read_id)
        if all(is_present(tags.get(tag)) for tag in ("L1", "L2", "L3")):
            self.stage_sets["ligation"].add(read_id)
        if is_present(tags.get("CB")):
            self.stage_sets["cb"].add(read_id)

        if self.modality == "rna":
            if is_present(tags.get("UM")):
                self.stage_sets["umi"].add(read_id)
        elif is_present(tags.get("MO")):
            self.stage_sets["modality"].add(read_id)

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
    parser.add_argument("--min-read-pairs-per-cell", default=100, type=int)
    return parser.parse_args()


def is_present(value: Optional[str]) -> bool:
    if value is None:
        return False
    cleaned = str(value).strip()
    return cleaned not in {"", "-", "NoMatch", "nomatch", "NA", "N/A", "None", "null"}


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


def add_warning(warnings: List[WarningRecord], unit: UnitReport, message: str) -> None:
    unit.add_warning(message)
    warnings.append(WarningRecord(unit.modality, unit.name, message))


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
    except Exception as exc:  # pragma: no cover - integration behavior
        add_warning(warnings, sample_unit, f"Could not parse tag records {path}: {exc}")


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


def get_target_units(
    units: Dict[Tuple[str, str], UnitReport],
    modality: str,
    split_name: str,
    suffix: str,
    sample_ids: Iterable[str],
) -> List[UnitReport]:
    sample, group, mark = parse_split_name(split_from_bam(Path(split_name), suffix), sample_ids, modality)
    target_units = [get_unit(units, modality, sample, "sample")]
    if group:
        target_units.append(get_unit(units, modality, f"{sample}_{group}", "sample_group"))
    if modality == "dna" and group and mark:
        target_units.append(get_unit(units, modality, f"{sample}_{group}_{mark}", "sample_group_mark"))
    return target_units


@dataclass
class BamObservation:
    mapped_ids: Set[str] = field(default_factory=set)
    gx_ids: Set[str] = field(default_factory=set)
    cb100_ids: Set[str] = field(default_factory=set)
    nonduplicate_ids: Set[str] = field(default_factory=set)
    cb100_available: bool = False
    warnings: List[str] = field(default_factory=list)


def read_cell_barcode(read) -> Tuple[Optional[str], Optional[str]]:
    if read.has_tag("CB"):
        value = str(read.get_tag("CB"))
        if is_present(value):
            return value, "CB"
    if read.has_tag("RG"):
        value = str(read.get_tag("RG"))
        if is_present(value):
            return value, "RG"
    return None, None


def collect_bam_observation(path: Path, min_read_pairs_per_cell: int) -> Tuple[Optional[BamObservation], Optional[str]]:
    try:
        import pysam  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on runtime env
        return None, f"pysam is unavailable; skipping BAM-derived categories for {path}: {exc}"

    observation = BamObservation()
    cb_to_read_names: Dict[str, Set[str]] = defaultdict(set)
    used_barcode_source: Optional[str] = None

    try:
        with pysam.AlignmentFile(str(path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                read_name = read.query_name
                observation.mapped_ids.add(read_name)

                if read.has_tag("GX") and is_present(str(read.get_tag("GX"))):
                    observation.gx_ids.add(read_name)

                if not read.is_duplicate:
                    observation.nonduplicate_ids.add(read_name)

                cell_barcode, barcode_source = read_cell_barcode(read)
                if cell_barcode:
                    cb_to_read_names[cell_barcode].add(read_name)
                    used_barcode_source = used_barcode_source or barcode_source
    except Exception as exc:
        return None, f"Could not read BAM-derived categories from {path}; skipping those categories: {exc}"

    if cb_to_read_names:
        observation.cb100_available = True
        high_count_barcodes = {
            cb for cb, read_names in cb_to_read_names.items() if len(read_names) >= min_read_pairs_per_cell
        }
        observation.cb100_ids = {
            read_name
            for cb in high_count_barcodes
            for read_name in cb_to_read_names[cb]
        }
        if used_barcode_source == "RG":
            observation.warnings.append(
                f"{path.name}: CB tag unavailable for at least one read; used RG tag as cell-barcode fallback for CB>100 +"
            )
    else:
        observation.warnings.append(
            f"{path.name}: no CB or RG tag found on mapped reads; omitting CB>100 + for affected units"
        )

    return observation, None


def add_rna_bam(
    units: Dict[Tuple[str, str], UnitReport],
    path: Path,
    sample_ids: Iterable[str],
    min_read_pairs_per_cell: int,
    warnings: List[WarningRecord],
) -> None:
    target_units = get_target_units(units, "rna", path.name, ".filtered_cells.bam", sample_ids)
    observation, warning = collect_bam_observation(path, min_read_pairs_per_cell)
    if warning:
        for unit in target_units:
            add_warning(warnings, unit, warning)
        return

    assert observation is not None
    for unit in target_units:
        unit.add_stage_ids("mapped", observation.mapped_ids)
        unit.add_stage_ids("gx", observation.gx_ids)
        if observation.cb100_available:
            unit.add_stage_ids("cb100", observation.cb100_ids)
        for message in observation.warnings:
            add_warning(warnings, unit, message)


def add_dna_markeddup_bam(
    units: Dict[Tuple[str, str], UnitReport],
    path: Path,
    sample_ids: Iterable[str],
    min_read_pairs_per_cell: int,
    warnings: List[WarningRecord],
) -> None:
    target_units = get_target_units(units, "dna", path.name, "_MarkedDup.bam", sample_ids)
    observation, warning = collect_bam_observation(path, min_read_pairs_per_cell)
    if warning:
        for unit in target_units:
            add_warning(warnings, unit, warning)
        return

    assert observation is not None
    for unit in target_units:
        unit.add_stage_ids("mapped", observation.mapped_ids)
        if observation.cb100_available:
            unit.add_stage_ids("cb100", observation.cb100_ids)
        unit.add_stage_ids("_unique_fallback", observation.nonduplicate_ids)
        for message in observation.warnings:
            add_warning(warnings, unit, message)


def add_dna_nodup_bam(
    units: Dict[Tuple[str, str], UnitReport],
    path: Path,
    sample_ids: Iterable[str],
    min_read_pairs_per_cell: int,
    warnings: List[WarningRecord],
) -> None:
    target_units = get_target_units(units, "dna", path.name, "_NoDup.bam", sample_ids)
    observation, warning = collect_bam_observation(path, min_read_pairs_per_cell)
    if warning:
        for unit in target_units:
            add_warning(warnings, unit, warning)
        return

    assert observation is not None
    for unit in target_units:
        unit.add_stage_ids("unique", observation.mapped_ids)


def apply_dna_unique_fallback(units: Dict[Tuple[str, str], UnitReport], warnings: List[WarningRecord]) -> None:
    for unit in units.values():
        if unit.modality != "dna":
            continue
        if "unique" in unit.stage_available or "_unique_fallback" not in unit.stage_available:
            continue
        unit.add_stage_ids("unique", set(unit.stage_sets.get("_unique_fallback", set())))
        add_warning(
            warnings,
            unit,
            "No readable DNA NoDup BAM was available for this unit; Unique + uses non-duplicate reads from MarkedDup BAM",
        )


def category_order(unit: UnitReport, min_read_pairs_per_cell: int) -> List[Tuple[str, str]]:
    categories = RNA_CATEGORIES if unit.modality == "rna" else DNA_CATEGORIES
    return [
        (stage, f"CB>{min_read_pairs_per_cell} +" if stage == "cb100" else label)
        for stage, label in categories
    ]


def write_placeholder_pdf(path: Path, title: str, message: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.axis("off")
        ax.text(0.02, 0.78, title, fontsize=16, weight="bold", transform=ax.transAxes)
        ax.text(0.02, 0.5, message, fontsize=11, transform=ax.transAxes, wrap=True)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)
    except Exception:
        path.write_text(f"{title}\n\n{message}\n", encoding="utf-8")


def write_upset_pdf(
    unit: UnitReport,
    pdf_path: Path,
    min_read_pairs_per_cell: int,
    warnings: List[WarningRecord],
) -> None:
    title = unit.name.replace("_", " ")
    categories = [
        (stage, label)
        for stage, label in category_order(unit, min_read_pairs_per_cell)
        if stage in unit.stage_available
    ]

    if not unit.tagged_ids:
        message = "No tag-record read identifiers were available for this unit."
        add_warning(warnings, unit, message)
        write_placeholder_pdf(pdf_path, title, message)
        return
    if len(categories) < 2:
        message = "Not enough available sequencing-efficiency categories to draw an UpSet plot."
        add_warning(warnings, unit, message)
        write_placeholder_pdf(pdf_path, title, message)
        return

    try:
        import pandas as pd
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from upsetplot import UpSet

        universe = sorted(unit.tagged_ids)
        labels = [label for _stage, label in categories]
        data = {
            label: [read_id in unit.stage_sets.get(stage, set()) for read_id in universe]
            for stage, label in categories
        }
        frame = pd.DataFrame(data)
        upset_data = frame.assign(read_records=1).groupby(labels, sort=False)["read_records"].sum()

        fig = plt.figure(figsize=(10, 12))
        try:
            upset = UpSet(
                upset_data,
                subset_size="count",
                sort_by="-degree",
                show_counts=True,
                show_percentages=True,
                min_subset_size="0.5%",
                facecolor="darkblue",
            )
        except TypeError:
            upset = UpSet(
                upset_data,
                subset_size="count",
                sort_by="-degree",
                show_counts=True,
                show_percentages=True,
                facecolor="darkblue",
            )
        upset.plot(fig=fig)
        fig.suptitle(title, fontsize=16, fontweight="bold", y=0.995)
        for axis in fig.axes:
            axis.tick_params(axis="both", labelsize=9)
            if axis.get_xlabel():
                axis.set_xlabel(axis.get_xlabel(), fontsize=10)
            if axis.get_ylabel():
                axis.set_ylabel(axis.get_ylabel(), fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
    except Exception as exc:
        message = f"UpSet plot could not be rendered for {unit.name}: {exc}"
        add_warning(warnings, unit, message)
        write_placeholder_pdf(pdf_path, title, message)


def emit_warnings(warnings: Sequence[WarningRecord]) -> None:
    seen: Set[WarningRecord] = set()
    for warning in warnings:
        if warning in seen:
            continue
        seen.add(warning)
        print(f"WARNING [{warning.modality}:{warning.unit}] {warning.message}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.min_read_pairs_per_cell < 1:
        print("--min-read-pairs-per-cell must be >= 1", file=sys.stderr)
        return 2

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
        add_rna_bam(units, path, sample_ids, args.min_read_pairs_per_cell, warnings)
    for path in args.dna_markeddup_bams:
        add_dna_markeddup_bam(units, path, sample_ids, args.min_read_pairs_per_cell, warnings)
    for path in args.dna_nodup_bams:
        add_dna_nodup_bam(units, path, sample_ids, args.min_read_pairs_per_cell, warnings)
    apply_dna_unique_fallback(units, warnings)

    for unit in sorted(units.values(), key=lambda item: (item.modality, item.name)):
        if not unit.tagged_ids:
            continue
        prefix = args.outdir / f"{unit.name}.{unit.modality}_sequencing_efficiency"
        write_upset_pdf(unit, Path(f"{prefix}.upset.pdf"), args.min_read_pairs_per_cell, warnings)
        for message in sorted(unit.warnings):
            warnings.append(WarningRecord(unit.modality, unit.name, message))

    emit_warnings(warnings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
