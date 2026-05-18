#!/usr/bin/env python3
"""Build RNA/DNA sequencing-efficiency UpSet PDFs from TrESFlow outputs.

The production path is an exact disk-backed reducer:
1. stream tag-record and BAM-derived category observations as integer bitmasks;
2. external-sort by unit/read id;
3. reduce sorted observations to per-combination counts;
4. render UpSet PDFs from aggregate counts.

This avoids keeping every read id in Python sets or materializing a per-read
pandas boolean matrix.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
import warnings as py_warnings
from collections import Counter, defaultdict
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

CATEGORY_TABLE = {"rna": RNA_CATEGORIES, "dna": DNA_CATEGORIES}
CATEGORY_BITS = {
    modality: {stage: 1 << idx for idx, (stage, _label) in enumerate(categories)}
    for modality, categories in CATEGORY_TABLE.items()
}
TAG_DERIVED_STAGES = {
    "rna": {"reads", "sample", "ligation", "cb", "umi"},
    "dna": {"reads", "sample", "ligation", "cb", "modality"},
}


@dataclass(frozen=True)
class WarningRecord:
    modality: str
    unit: str
    message: str


@dataclass
class UnitState:
    name: str
    modality: str
    level: str
    stage_available: Set[str] = field(default_factory=set)
    warnings: Set[str] = field(default_factory=set)
    combination_counts: Counter[int] = field(default_factory=Counter)
    fallback_unique_observed: bool = False
    unique_bam_available: bool = False

    def add_warning(self, message: str) -> None:
        self.warnings.add(message)

    def register(self, *stages: str) -> None:
        self.stage_available.update(stage for stage in stages if stage in CATEGORY_BITS[self.modality])


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
    parser.add_argument("--sort-parallel", default=1, type=int)
    parser.add_argument("--sort-buffer", default="2G")
    parser.add_argument("--tmpdir", type=Path)
    parser.add_argument("--debug-counts", type=Path, help=argparse.SUPPRESS)
    return parser.parse_args()


def log(message: str) -> None:
    print(f"[sequencing-efficiency] {message}", file=sys.stderr, flush=True)


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


def get_unit(units: Dict[Tuple[str, str], UnitState], modality: str, name: str, level: str) -> UnitState:
    key = unit_key(modality, name)
    if key not in units:
        units[key] = UnitState(name=name, modality=modality, level=level)
    return units[key]


def add_warning(warnings: List[WarningRecord], unit: UnitState, message: str) -> None:
    unit.add_warning(message)
    warnings.append(WarningRecord(unit.modality, unit.name, message))


def category_bit(modality: str, stage: str) -> int:
    return CATEGORY_BITS[modality][stage]


def tag_mask(modality: str, tags: Dict[str, str]) -> int:
    mask = category_bit(modality, "reads")
    if is_present(tags.get("SB")):
        mask |= category_bit(modality, "sample")
    if all(is_present(tags.get(tag)) for tag in ("L1", "L2", "L3")):
        mask |= category_bit(modality, "ligation")
    if is_present(tags.get("CB")):
        mask |= category_bit(modality, "cb")
    if modality == "rna":
        if is_present(tags.get("UM")):
            mask |= category_bit(modality, "umi")
    elif is_present(tags.get("MO")):
        mask |= category_bit(modality, "modality")
    return mask


def write_observation(handle, unit: UnitState, read_id: str, mask: int) -> None:
    if "\t" in read_id:
        read_id = read_id.replace("\t", " ")
    handle.write(f"{unit.modality}\t{unit.name}\t{read_id}\t{mask}\n")


def add_tag_records(
    units: Dict[Tuple[str, str], UnitState],
    modality: str,
    path: Path,
    sb_maps: Dict[str, Dict[str, str]],
    mo_maps: Dict[str, Dict[Tuple[str, str], str]],
    warnings: List[WarningRecord],
    obs_handle,
) -> Tuple[int, int]:
    sample = sample_from_tag_record(path, modality)
    sample_unit = get_unit(units, modality, sample, "sample")
    sample_unit.register(*TAG_DERIVED_STAGES[modality])
    started = time.monotonic()
    line_count = 0
    obs_count = 0

    try:
        with open_maybe_gzip(path) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                read_id, tags = parse_tag_record_line(raw_line)
                if read_id is None:
                    continue

                mask = tag_mask(modality, tags)
                write_observation(obs_handle, sample_unit, read_id, mask)
                obs_count += 1

                group = resolve_group(sample, tags.get("SB"), sb_maps)
                if group:
                    group_unit = get_unit(units, modality, f"{sample}_{group}", "sample_group")
                    group_unit.register(*TAG_DERIVED_STAGES[modality])
                    write_observation(obs_handle, group_unit, read_id, mask)
                    obs_count += 1
                    if modality == "dna":
                        mark = resolve_mark(sample, group, tags.get("MO"), mo_maps)
                        if mark:
                            mark_unit = get_unit(units, modality, f"{sample}_{group}_{mark}", "sample_group_mark")
                            mark_unit.register(*TAG_DERIVED_STAGES[modality])
                            write_observation(obs_handle, mark_unit, read_id, mask)
                            obs_count += 1

                line_count += 1
                if line_count % 5_000_000 == 0:
                    elapsed = time.monotonic() - started
                    log(f"{path.name}: streamed {line_count:,} tag records in {elapsed:.1f}s")
    except Exception as exc:  # pragma: no cover - integration behavior
        add_warning(warnings, sample_unit, f"Could not parse tag records {path}: {exc}")

    elapsed = time.monotonic() - started
    log(f"{path.name}: wrote {obs_count:,} exact tag observations from {line_count:,} records in {elapsed:.1f}s")
    return line_count, obs_count


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
    units: Dict[Tuple[str, str], UnitState],
    modality: str,
    split_name: str,
    suffix: str,
    sample_ids: Iterable[str],
) -> List[UnitState]:
    sample, group, mark = parse_split_name(split_from_bam(Path(split_name), suffix), sample_ids, modality)
    target_units = [get_unit(units, modality, sample, "sample")]
    if group:
        target_units.append(get_unit(units, modality, f"{sample}_{group}", "sample_group"))
    if modality == "dna" and group and mark:
        target_units.append(get_unit(units, modality, f"{sample}_{group}_{mark}", "sample_group_mark"))
    return target_units


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


def run_sort(
    input_path: Path,
    output_path: Path,
    tmpdir: Path,
    parallel: int,
    buffer_size: str,
    key_args: Sequence[str],
    unique: bool = False,
) -> None:
    cmd = [
        "sort",
        "-t",
        "\t",
        "--parallel",
        str(max(1, parallel)),
        "-S",
        buffer_size,
        "-T",
        str(tmpdir),
    ]
    if unique:
        cmd.append("-u")
    cmd.extend(key_args)
    cmd.extend([str(input_path), "-o", str(output_path)])
    started = time.monotonic()
    subprocess.run(cmd, check=True)
    log(f"sorted {input_path.name} in {time.monotonic() - started:.1f}s")


def emit_cb100_observations(
    sorted_pairs_path: Path,
    target_units: Sequence[UnitState],
    modality: str,
    threshold: int,
    obs_handle,
) -> int:
    cb100_bit = category_bit(modality, "cb100")
    emitted = 0
    current_cb: Optional[str] = None
    read_names: List[str] = []

    def flush_current() -> int:
        if current_cb is None or len(read_names) < threshold:
            return 0
        local_count = 0
        for read_name in read_names:
            for unit in target_units:
                write_observation(obs_handle, unit, read_name, cb100_bit)
                local_count += 1
        return local_count

    with open(sorted_pairs_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            cb, read_name = raw_line.rstrip("\n").split("\t", 1)
            if cb != current_cb:
                emitted += flush_current()
                current_cb = cb
                read_names = [read_name]
            else:
                read_names.append(read_name)
    emitted += flush_current()
    return emitted


def process_bam(
    units: Dict[Tuple[str, str], UnitState],
    path: Path,
    modality: str,
    bam_kind: str,
    suffix: str,
    sample_ids: Iterable[str],
    min_read_pairs_per_cell: int,
    warnings: List[WarningRecord],
    obs_handle,
    fallback_handle,
    tmpdir: Path,
    sort_parallel: int,
    sort_buffer: str,
) -> None:
    target_units = get_target_units(units, modality, path.name, suffix, sample_ids)
    started = time.monotonic()

    try:
        import pysam  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on runtime env
        for unit in target_units:
            add_warning(warnings, unit, f"pysam is unavailable; skipping BAM-derived categories for {path}: {exc}")
        return

    if bam_kind == "rna_filtered":
        bam_stages = ("mapped", "gx")
        mapped_bit = category_bit("rna", "mapped")
        gx_bit = category_bit("rna", "gx")
        cb100_pair_path = tmpdir / f"{path.name}.cb_pairs.tsv"
    elif bam_kind == "dna_markeddup":
        bam_stages = ("mapped",)
        mapped_bit = category_bit("dna", "mapped")
        unique_bit = category_bit("dna", "unique")
        cb100_pair_path = tmpdir / f"{path.name}.cb_pairs.tsv"
    elif bam_kind == "dna_nodup":
        bam_stages = ("unique",)
        unique_bit = category_bit("dna", "unique")
        cb100_pair_path = None
    else:  # pragma: no cover - guarded by callers
        raise ValueError(f"Unsupported BAM kind: {bam_kind}")

    reads_seen = 0
    mapped_seen = 0
    category_observations = 0
    cb_pair_count = 0
    used_rg = False
    cb_pair_handle = None

    try:
        if cb100_pair_path is not None:
            cb_pair_handle = open(cb100_pair_path, "wt", encoding="utf-8")
        with pysam.AlignmentFile(str(path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                reads_seen += 1
                if read.is_unmapped:
                    continue
                mapped_seen += 1
                read_name = read.query_name

                if bam_kind in {"rna_filtered", "dna_markeddup"}:
                    for unit in target_units:
                        write_observation(obs_handle, unit, read_name, mapped_bit)
                        category_observations += 1

                if bam_kind == "rna_filtered" and read.has_tag("GX") and is_present(str(read.get_tag("GX"))):
                    for unit in target_units:
                        write_observation(obs_handle, unit, read_name, gx_bit)
                        category_observations += 1

                if bam_kind == "dna_markeddup" and not read.is_duplicate:
                    for unit in target_units:
                        write_observation(fallback_handle, unit, read_name, unique_bit)
                        unit.fallback_unique_observed = True

                if bam_kind == "dna_nodup":
                    for unit in target_units:
                        write_observation(obs_handle, unit, read_name, unique_bit)
                        category_observations += 1

                if cb_pair_handle is not None:
                    cell_barcode, barcode_source = read_cell_barcode(read)
                    if cell_barcode:
                        if barcode_source == "RG":
                            used_rg = True
                        cb_pair_handle.write(f"{cell_barcode}\t{read_name}\n")
                        cb_pair_count += 1

                if reads_seen % 5_000_000 == 0:
                    log(f"{path.name}: streamed {reads_seen:,} BAM records in {time.monotonic() - started:.1f}s")
    except Exception as exc:
        for unit in target_units:
            add_warning(warnings, unit, f"Could not read BAM-derived categories from {path}; skipping those categories: {exc}")
        return
    finally:
        if cb_pair_handle is not None:
            cb_pair_handle.close()

    for unit in target_units:
        unit.register(*bam_stages)
        if bam_kind == "dna_nodup":
            unit.unique_bam_available = True

    cb100_observations = 0
    if cb100_pair_path is not None:
        if cb_pair_count:
            sorted_pairs_path = tmpdir / f"{path.name}.cb_pairs.sorted.tsv"
            run_sort(
                cb100_pair_path,
                sorted_pairs_path,
                tmpdir,
                sort_parallel,
                sort_buffer,
                ["-k1,1", "-k2,2"],
                unique=True,
            )
            cb100_observations = emit_cb100_observations(
                sorted_pairs_path,
                target_units,
                modality,
                min_read_pairs_per_cell,
                obs_handle,
            )
            for unit in target_units:
                unit.register("cb100")
            if used_rg:
                for unit in target_units:
                    add_warning(
                        warnings,
                        unit,
                        f"{path.name}: CB tag unavailable for at least one read; used RG tag as cell-barcode fallback for CB>{min_read_pairs_per_cell} +",
                    )
            cb100_pair_path.unlink(missing_ok=True)
            sorted_pairs_path.unlink(missing_ok=True)
        else:
            cb100_pair_path.unlink(missing_ok=True)
            for unit in target_units:
                add_warning(
                    warnings,
                    unit,
                    f"{path.name}: no CB or RG tag found on mapped reads; omitting CB>{min_read_pairs_per_cell} + for affected units",
                )

    elapsed = time.monotonic() - started
    log(
        f"{path.name}: scanned {reads_seen:,} BAM records, {mapped_seen:,} mapped records, "
        f"wrote {category_observations + cb100_observations:,} observations in {elapsed:.1f}s"
    )


def append_needed_unique_fallbacks(
    fallback_path: Path,
    obs_handle,
    units: Dict[Tuple[str, str], UnitState],
    warnings: List[WarningRecord],
) -> int:
    needed = {
        (unit.modality, unit.name)
        for unit in units.values()
        if unit.modality == "dna" and not unit.unique_bam_available and unit.fallback_unique_observed
    }
    if not needed or not fallback_path.exists() or fallback_path.stat().st_size == 0:
        return 0

    appended = 0
    with open(fallback_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            modality, unit_name, _read_id, _mask = raw_line.rstrip("\n").split("\t", 3)
            if (modality, unit_name) not in needed:
                continue
            obs_handle.write(raw_line)
            appended += 1

    for modality, unit_name in sorted(needed):
        unit = units[(modality, unit_name)]
        unit.register("unique")
        add_warning(
            warnings,
            unit,
            "No readable DNA NoDup BAM was available for this unit; Unique + uses non-duplicate reads from MarkedDup BAM",
        )
    log(f"appended {appended:,} DNA Unique + fallback observations from MarkedDup BAMs")
    return appended


def reduce_sorted_observations(
    sorted_obs_path: Path,
    units: Dict[Tuple[str, str], UnitState],
) -> int:
    started = time.monotonic()
    reduced_reads = 0
    current_key: Optional[Tuple[str, str, str]] = None
    current_mask = 0

    def flush_current() -> None:
        nonlocal reduced_reads, current_mask, current_key
        if current_key is None:
            return
        modality, unit_name, _read_id = current_key
        reads_bit = category_bit(modality, "reads")
        if current_mask & reads_bit:
            unit = get_unit(units, modality, unit_name, "unknown")
            unit.combination_counts[current_mask] += 1
            reduced_reads += 1

    with open(sorted_obs_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            modality, unit_name, read_id, mask_text = raw_line.rstrip("\n").split("\t", 3)
            key = (modality, unit_name, read_id)
            mask = int(mask_text)
            if key != current_key:
                flush_current()
                current_key = key
                current_mask = mask
            else:
                current_mask |= mask
    flush_current()

    log(f"reduced {reduced_reads:,} tag-record read ids in {time.monotonic() - started:.1f}s")
    return reduced_reads


def category_order(unit: UnitState, min_read_pairs_per_cell: int) -> List[Tuple[str, str]]:
    return [
        (stage, f"CB>{min_read_pairs_per_cell} +" if stage == "cb100" else label)
        for stage, label in CATEGORY_TABLE[unit.modality]
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
    unit: UnitState,
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

    if not unit.combination_counts:
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

        labels = [label for _stage, label in categories]
        bit_order = [category_bit(unit.modality, stage) for stage, _label in categories]
        collapsed: Counter[Tuple[bool, ...]] = Counter()
        for mask, count in unit.combination_counts.items():
            collapsed[tuple(bool(mask & bit) for bit in bit_order)] += count

        index = pd.MultiIndex.from_tuples(list(collapsed.keys()), names=labels)
        upset_data = pd.Series(list(collapsed.values()), index=index, name="read_records")

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
        with py_warnings.catch_warnings():
            py_warnings.filterwarnings("ignore", message="This figure includes Axes that are not compatible with tight_layout")
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


def write_debug_counts(path: Path, units: Dict[Tuple[str, str], UnitState]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt", encoding="utf-8") as handle:
        handle.write("modality\tunit\tmask\tread_records\tcategories\n")
        for unit in sorted(units.values(), key=lambda item: (item.modality, item.name)):
            for mask, count in sorted(unit.combination_counts.items()):
                categories = [
                    stage
                    for stage, _label in CATEGORY_TABLE[unit.modality]
                    if mask & category_bit(unit.modality, stage)
                ]
                handle.write(f"{unit.modality}\t{unit.name}\t{mask}\t{count}\t{','.join(categories)}\n")


def prepare_temp_dir(args: argparse.Namespace) -> Path:
    root = args.tmpdir if args.tmpdir else args.outdir
    root.mkdir(parents=True, exist_ok=True)
    tmpdir = root / f".sequencing_efficiency_tmp_{time.strftime('%Y%m%d%H%M%S')}_{os.getpid()}"
    tmpdir.mkdir(parents=True, exist_ok=False)
    return tmpdir


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.min_read_pairs_per_cell < 1:
        print("--min-read-pairs-per-cell must be >= 1", file=sys.stderr)
        return 2
    if args.sort_parallel < 1:
        print("--sort-parallel must be >= 1", file=sys.stderr)
        return 2

    rna_sb_paths, dna_sb_paths = split_sb_group_map_paths(args.sb_group_maps)
    rna_sb_maps, rna_groups_by_sample = load_sb_group_maps(rna_sb_paths)
    dna_sb_maps, dna_groups_by_sample = load_sb_group_maps(dna_sb_paths)
    mo_maps = load_dna_mo_maps(args.dna_mo_maps)
    units: Dict[Tuple[str, str], UnitState] = {}
    warnings: List[WarningRecord] = []
    tmpdir = prepare_temp_dir(args)
    obs_path = tmpdir / "observations.tsv"
    fallback_path = tmpdir / "dna_unique_fallback.tsv"
    sorted_obs_path = tmpdir / "observations.sorted.tsv"

    try:
        for sample, groups in rna_groups_by_sample.items():
            for group in groups:
                get_unit(units, "rna", f"{sample}_{group}", "sample_group")
        for sample, groups in dna_groups_by_sample.items():
            for group in groups:
                get_unit(units, "dna", f"{sample}_{group}", "sample_group")

        with open(obs_path, "wt", encoding="utf-8") as obs_handle, open(
            fallback_path, "wt", encoding="utf-8"
        ) as fallback_handle:
            for path in args.rna_tag_records:
                add_tag_records(units, "rna", path, rna_sb_maps, mo_maps, warnings, obs_handle)
            for path in args.dna_tag_records:
                add_tag_records(units, "dna", path, dna_sb_maps, mo_maps, warnings, obs_handle)

            sample_ids = {sample_from_tag_record(path, "rna") for path in args.rna_tag_records}
            sample_ids.update(sample_from_tag_record(path, "dna") for path in args.dna_tag_records)
            sample_ids.update(rna_sb_maps.keys())
            sample_ids.update(dna_sb_maps.keys())

            for path in args.rna_filtered_bams:
                process_bam(
                    units,
                    path,
                    "rna",
                    "rna_filtered",
                    ".filtered_cells.bam",
                    sample_ids,
                    args.min_read_pairs_per_cell,
                    warnings,
                    obs_handle,
                    fallback_handle,
                    tmpdir,
                    args.sort_parallel,
                    args.sort_buffer,
                )
            for path in args.dna_markeddup_bams:
                process_bam(
                    units,
                    path,
                    "dna",
                    "dna_markeddup",
                    "_MarkedDup.bam",
                    sample_ids,
                    args.min_read_pairs_per_cell,
                    warnings,
                    obs_handle,
                    fallback_handle,
                    tmpdir,
                    args.sort_parallel,
                    args.sort_buffer,
                )
            for path in args.dna_nodup_bams:
                process_bam(
                    units,
                    path,
                    "dna",
                    "dna_nodup",
                    "_NoDup.bam",
                    sample_ids,
                    args.min_read_pairs_per_cell,
                    warnings,
                    obs_handle,
                    fallback_handle,
                    tmpdir,
                    args.sort_parallel,
                    args.sort_buffer,
                )

            append_needed_unique_fallbacks(fallback_path, obs_handle, units, warnings)

        if obs_path.stat().st_size > 0:
            run_sort(
                obs_path,
                sorted_obs_path,
                tmpdir,
                args.sort_parallel,
                args.sort_buffer,
                ["-k1,1", "-k2,2", "-k3,3"],
                unique=False,
            )
            reduce_sorted_observations(sorted_obs_path, units)
        else:
            log("no sequencing-efficiency observations were emitted")

        for unit in sorted(units.values(), key=lambda item: (item.modality, item.name)):
            if not unit.combination_counts:
                continue
            prefix = args.outdir / f"{unit.name}.{unit.modality}_sequencing_efficiency"
            write_upset_pdf(unit, Path(f"{prefix}.upset.pdf"), args.min_read_pairs_per_cell, warnings)
            for message in sorted(unit.warnings):
                warnings.append(WarningRecord(unit.modality, unit.name, message))

        if args.debug_counts:
            write_debug_counts(args.debug_counts, units)

        emit_warnings(warnings)
        return 0
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
