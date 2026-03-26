#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd


RNA_STAGE_ORDER = [
    "start_reads",
    "all_proper_barcodes_reads",
    "usable_in_genes_reads",
    "unique_in_genes_reads",
]

DNA_STAGE_ORDER = [
    "start_reads",
    "all_proper_barcodes_reads",
    "properly_mapped_reads",
    "unique_nodup_reads",
]

DNA_MARK_STAGE_ORDER = [
    "all_proper_barcodes_reads",
    "properly_mapped_reads",
    "unique_nodup_reads",
]

STAGE_LABELS = {
    "start_reads": "Start",
    "all_proper_barcodes_reads": "All Proper Barcodes",
    "usable_in_genes_reads": "Usable In Genes",
    "unique_in_genes_reads": "Unique In Genes",
    "properly_mapped_reads": "Properly Mapped",
    "unique_nodup_reads": "Unique NoDup",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate TrESFlow core QC summaries and plots.")
    subparsers = parser.add_subparsers(dest="mode", required=True)

    rna = subparsers.add_parser("rna", help="Build RNA QC outputs")
    rna.add_argument("--sb-group-map", required=True)
    rna.add_argument("--sample-counts", action="append", default=[])
    rna.add_argument("--cell-counts", action="append", default=[])
    rna.add_argument("--solo-dir", action="append", default=[])
    rna.add_argument("--outdir", required=True)

    dna = subparsers.add_parser("dna", help="Build DNA QC outputs")
    dna.add_argument("--sb-group-map", required=True)
    dna.add_argument("--dna-mo-map", required=True)
    dna.add_argument("--sample-counts", action="append", default=[])
    dna.add_argument("--tag-records", action="append", default=[])
    dna.add_argument("--aligned-bam", action="append", default=[])
    dna.add_argument("--nodup-bam", action="append", default=[])
    dna.add_argument("--outdir", required=True)

    return parser.parse_args()


def load_sb_group_map(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]], dict[str, dict[str, list[str]]]]:
    df = pd.read_csv(path, sep="\t")
    required = {"sample", "sb_group", "sb_bc"}
    if not required.issubset(df.columns):
        raise ValueError(f"Invalid sb_group_map columns in {path}: {sorted(df.columns)}")

    sample_order = list(dict.fromkeys(df["sample"].tolist()))
    group_order = {}
    barcodes = {}
    for sample in sample_order:
        sample_df = df[df["sample"] == sample]
        groups = list(dict.fromkeys(sample_df["sb_group"].tolist()))
        group_order[sample] = groups
        barcodes[sample] = {
            group: sample_df[sample_df["sb_group"] == group]["sb_bc"].tolist()
            for group in groups
        }

    return df, group_order, barcodes


def load_dna_mo_map(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]], dict[tuple[str, str, str], str]]:
    df = pd.read_csv(path, sep="\t")
    required = {"sample", "sb_group", "mark", "mo_bc"}
    if not required.issubset(df.columns):
        raise ValueError(f"Invalid dna_mo_map columns in {path}: {sorted(df.columns)}")

    mark_order = {}
    mark_lookup = {}
    for sample in dict.fromkeys(df["sample"].tolist()):
        sample_df = df[df["sample"] == sample]
        mark_order[sample] = list(dict.fromkeys(sample_df["mark"].tolist()))
        for row in sample_df.itertuples(index=False):
            mark_lookup[(row.sample, row.sb_group, row.mo_bc)] = row.mark

    return df, mark_order, mark_lookup


def read_counts_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=["reads", "barcode"])
    df["reads"] = df["reads"].astype(int)
    df["barcode"] = df["barcode"].astype(str)
    return df


def read_tag_records(path: Path) -> pd.DataFrame:
    rows = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            fields = {}
            for token in parts[1:]:
                if token.startswith("SB:Z:"):
                    fields["sb"] = token[len("SB:Z:") :]
                elif token.startswith("MO:Z:"):
                    fields["mo"] = token[len("MO:Z:") :]
                elif token.startswith("CB:Z:"):
                    fields["cb"] = token[len("CB:Z:") :]
            rows.append(fields)

    df = pd.DataFrame(rows)
    for column in ("sb", "mo", "cb"):
        if column not in df.columns:
            df[column] = ""
        df[column] = df[column].fillna("").astype(str)
    return df


def extract_sample_from_filename(path: Path, suffix: str) -> str:
    name = path.name
    if not name.endswith(suffix):
        raise ValueError(f"Unexpected file name for suffix '{suffix}': {path}")
    return name[: -len(suffix)]


def match_sample_prefix(name: str, sample_ids: list[str]) -> tuple[str, str]:
    for sample in sorted(sample_ids, key=len, reverse=True):
        prefix = f"{sample}_"
        if name.startswith(prefix):
            return sample, name[len(prefix):]
    raise ValueError(f"Could not match sample prefix for '{name}'")


def parse_rna_split_name(path: Path, sample_ids: list[str]) -> tuple[str, str]:
    stem = path.name
    if stem.endswith(".Solo.outGeneFull"):
        stem = stem[: -len(".Solo.outGeneFull")]
    sample, group = match_sample_prefix(stem, sample_ids)
    return sample, group


def parse_dna_bam_name(path: Path, sample_ids: list[str], mark_order: dict[str, list[str]], suffix: str) -> tuple[str, str, str]:
    stem = path.name
    if not stem.endswith(suffix):
        raise ValueError(f"Unexpected DNA BAM name '{path.name}' for suffix '{suffix}'")
    stem = stem[: -len(suffix)]
    sample, remainder = match_sample_prefix(stem, sample_ids)
    for mark in sorted(mark_order.get(sample, []), key=len, reverse=True):
        mark_suffix = f"_{mark}"
        if remainder.endswith(mark_suffix):
            group = remainder[: -len(mark_suffix)]
            if not group:
                break
            return sample, group, mark
    raise ValueError(f"Could not derive DNA group and mark from '{path.name}'")


def integer_formatter() -> FuncFormatter:
    return FuncFormatter(lambda value, _pos: f"{int(value):,}")


def annotate_bars(ax: plt.Axes) -> None:
    for container in ax.containers:
        labels = []
        for bar in container:
            height = bar.get_height()
            labels.append("" if height <= 0 else f"{int(height):,}")
        ax.bar_label(container, labels=labels, fontsize=8, padding=2, rotation=90)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_table(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        raise ValueError(f"Refusing to write empty QC table: {path}")
    df.to_csv(path, sep="\t", index=False)


def plot_group_overview(df: pd.DataFrame, sample: str, group_col: str, stage_order: list[str], out_path: Path, title: str) -> None:
    plot_df = df[df["sample"] == sample].copy()
    plot_df["stage"] = pd.Categorical(plot_df["stage"], categories=stage_order, ordered=True)
    plot_df = plot_df.sort_values(["stage", group_col])

    pivot = plot_df.pivot(index="stage", columns=group_col, values="reads").fillna(0)
    pivot.index = [STAGE_LABELS[stage] for stage in pivot.index]

    fig, ax = plt.subplots(figsize=(10, 6))
    pivot.plot(kind="bar", ax=ax)
    ax.set_title(title)
    ax.set_ylabel("Reads")
    ax.set_xlabel("")
    ax.yaxis.set_major_formatter(integer_formatter())
    ax.legend(title=group_col.replace("_", " ").title())
    annotate_bars(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_single_group(df: pd.DataFrame, sample: str, group: str, stage_order: list[str], out_path: Path, title: str) -> None:
    plot_df = df[(df["sample"] == sample) & (df["group"] == group)].copy()
    plot_df["stage"] = pd.Categorical(plot_df["stage"], categories=stage_order, ordered=True)
    plot_df = plot_df.sort_values("stage")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar([STAGE_LABELS[stage] for stage in plot_df["stage"]], plot_df["reads"], color="#4c78a8")
    ax.set_title(title)
    ax.set_ylabel("Reads")
    ax.yaxis.set_major_formatter(integer_formatter())
    annotate_bars(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_group_marks(df: pd.DataFrame, sample: str, group: str, stage_order: list[str], out_path: Path, title: str) -> None:
    plot_df = df[(df["sample"] == sample) & (df["group"] == group)].copy()
    plot_df["stage"] = pd.Categorical(plot_df["stage"], categories=stage_order, ordered=True)
    plot_df = plot_df.sort_values(["stage", "mark"])

    pivot = plot_df.pivot(index="stage", columns="mark", values="reads").fillna(0)
    pivot.index = [STAGE_LABELS[stage] for stage in pivot.index]

    fig, ax = plt.subplots(figsize=(10, 6))
    pivot.plot(kind="bar", ax=ax)
    ax.set_title(title)
    ax.set_ylabel("Reads")
    ax.set_xlabel("")
    ax.yaxis.set_major_formatter(integer_formatter())
    ax.legend(title="Mark")
    annotate_bars(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def build_rna_qc(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    rna_plot_dir = outdir / "rna"
    group_plot_dir = rna_plot_dir / "groups"
    ensure_dir(rna_plot_dir)
    ensure_dir(group_plot_dir)

    _, group_order, group_barcodes = load_sb_group_map(Path(args.sb_group_map))
    sample_ids = list(group_order.keys())

    sample_counts_files = {
        extract_sample_from_filename(Path(path), ".sample_barcode.counts.tsv"): Path(path)
        for path in args.sample_counts
    }
    cell_counts_files = {
        extract_sample_from_filename(Path(path), ".cell.counts.tsv"): Path(path)
        for path in args.cell_counts
    }

    solo_dirs_by_sample = {}
    for raw_path in args.solo_dir:
        solo_dir = Path(raw_path)
        sample, group = parse_rna_split_name(solo_dir, sample_ids)
        solo_dirs_by_sample.setdefault(sample, {})[group] = solo_dir

    group_rows = []
    sample_rows = []

    for sample in sample_ids:
        if sample not in sample_counts_files:
            raise ValueError(f"Missing RNA sample counts for sample '{sample}'")
        if sample not in cell_counts_files:
            raise ValueError(f"Missing RNA cell counts for sample '{sample}'")

        sample_counts = read_counts_table(sample_counts_files[sample])
        cell_counts = read_counts_table(cell_counts_files[sample])
        sb_lookup = {}
        sb_length = None
        for group, barcodes in group_barcodes[sample].items():
            for barcode in barcodes:
                sb_lookup[barcode] = group
                sb_length = len(barcode) if sb_length is None else sb_length

        if sb_length is None:
            raise ValueError(f"No RNA sample barcodes defined for sample '{sample}'")

        group_stage_totals = {group: {stage: 0 for stage in RNA_STAGE_ORDER} for group in group_order[sample]}

        for group, barcodes in group_barcodes[sample].items():
            group_stage_totals[group]["start_reads"] = int(
                sample_counts[sample_counts["barcode"].isin(barcodes)]["reads"].sum()
            )

        for row in cell_counts.itertuples(index=False):
            prefix = row.barcode[:sb_length]
            if prefix not in sb_lookup:
                continue
            group_stage_totals[sb_lookup[prefix]]["all_proper_barcodes_reads"] += int(row.reads)

        sample_solo_dirs = solo_dirs_by_sample.get(sample, {})
        for group in group_order[sample]:
            if group not in sample_solo_dirs:
                raise ValueError(f"Missing RNA STARsolo directory for sample '{sample}' group '{group}'")

            cell_reads_stats = sample_solo_dirs[group] / "CellReads.stats"
            if not cell_reads_stats.exists():
                raise ValueError(f"Missing RNA CellReads.stats: {cell_reads_stats}")

            stats_df = pd.read_csv(cell_reads_stats, sep="\t")
            stats_df = stats_df[stats_df["CB"] != "CBnotInPasslist"].copy()
            for column in ("countedU", "countedM"):
                stats_df[column] = pd.to_numeric(stats_df[column], errors="raise")

            group_stage_totals[group]["usable_in_genes_reads"] = int((stats_df["countedU"] + stats_df["countedM"]).sum())
            group_stage_totals[group]["unique_in_genes_reads"] = int(stats_df["countedU"].sum())

        for group in group_order[sample]:
            for stage in RNA_STAGE_ORDER:
                group_rows.append(
                    {
                        "sample": sample,
                        "group": group,
                        "stage": stage,
                        "reads": int(group_stage_totals[group][stage]),
                    }
                )

        for stage in RNA_STAGE_ORDER:
            sample_rows.append(
                {
                    "sample": sample,
                    "stage": stage,
                    "reads": int(sum(group_stage_totals[group][stage] for group in group_order[sample])),
                }
            )

    group_df = pd.DataFrame(group_rows)
    sample_df = pd.DataFrame(sample_rows)
    write_table(group_df, outdir / "rna_group_stage_counts.tsv")
    write_table(sample_df, outdir / "rna_sample_stage_counts.tsv")

    for sample in sample_ids:
        plot_group_overview(
            group_df,
            sample,
            "group",
            RNA_STAGE_ORDER,
            rna_plot_dir / f"{sample}.group_stage_counts.png",
            f"RNA QC: {sample}",
        )
        for group in group_order[sample]:
            plot_single_group(
                group_df,
                sample,
                group,
                RNA_STAGE_ORDER,
                group_plot_dir / f"{sample}__{group}.png",
                f"RNA QC: {sample} / {group}",
            )

    return 0


def samtools_count(path: Path) -> int:
    samtools_bin = os.environ.get("SAMTOOLS_BIN", "samtools")
    result = subprocess.run(
        [samtools_bin, "view", "-c", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(result.stdout.strip())


def build_dna_qc(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    dna_plot_dir = outdir / "dna"
    group_plot_dir = dna_plot_dir / "groups"
    group_mark_plot_dir = dna_plot_dir / "group_marks"
    ensure_dir(dna_plot_dir)
    ensure_dir(group_plot_dir)
    ensure_dir(group_mark_plot_dir)

    _, group_order, group_barcodes = load_sb_group_map(Path(args.sb_group_map))
    _, mark_order, mark_lookup = load_dna_mo_map(Path(args.dna_mo_map))
    sample_ids = list(group_order.keys())

    sample_counts_files = {
        extract_sample_from_filename(Path(path), ".dna_sample_barcode.counts.tsv"): Path(path)
        for path in args.sample_counts
    }
    tag_record_files = {
        extract_sample_from_filename(Path(path), ".dna_tag_records.tsv"): Path(path)
        for path in args.tag_records
    }

    aligned_by_sample = {}
    for raw_path in args.aligned_bam:
        bam_path = Path(raw_path)
        sample, group, mark = parse_dna_bam_name(bam_path, sample_ids, mark_order, ".bam")
        aligned_by_sample.setdefault(sample, {})[(group, mark)] = bam_path

    nodup_by_sample = {}
    for raw_path in args.nodup_bam:
        bam_path = Path(raw_path)
        sample, group, mark = parse_dna_bam_name(bam_path, sample_ids, mark_order, "_NoDup.bam")
        nodup_by_sample.setdefault(sample, {})[(group, mark)] = bam_path

    group_rows = []
    group_mark_rows = []
    sample_rows = []

    for sample in sample_ids:
        if sample not in sample_counts_files:
            raise ValueError(f"Missing DNA sample counts for sample '{sample}'")
        if sample not in tag_record_files:
            raise ValueError(f"Missing DNA tag records for sample '{sample}'")

        sample_counts = read_counts_table(sample_counts_files[sample])
        tag_records = read_tag_records(tag_record_files[sample])

        sb_lookup = {}
        sb_length = None
        for group, barcodes in group_barcodes[sample].items():
            for barcode in barcodes:
                sb_lookup[barcode] = group
                sb_length = len(barcode) if sb_length is None else sb_length

        if sb_length is None:
            raise ValueError(f"No DNA sample barcodes defined for sample '{sample}'")

        sample_mark_lookup = {
            (group, mo_bc): mark
            for (lookup_sample, group, mo_bc), mark in mark_lookup.items()
            if lookup_sample == sample
        }
        group_stage_totals = {group: {stage: 0 for stage in DNA_STAGE_ORDER} for group in group_order[sample]}
        group_mark_stage_totals = {
            (group, mark): {stage: 0 for stage in DNA_MARK_STAGE_ORDER}
            for group in group_order[sample]
            for mark in mark_order.get(sample, [])
        }

        for group, barcodes in group_barcodes[sample].items():
            group_stage_totals[group]["start_reads"] = int(
                sample_counts[sample_counts["barcode"].isin(barcodes)]["reads"].sum()
            )

        for row in tag_records.itertuples(index=False):
            if row.sb not in sb_lookup:
                continue
            if row.mo == "NoMatch" or row.cb == "NoMatch":
                continue
            group = sb_lookup[row.sb]
            mark = sample_mark_lookup.get((group, row.mo))
            if mark is None:
                continue
            group_stage_totals[group]["all_proper_barcodes_reads"] += 1
            group_mark_stage_totals[(group, mark)]["all_proper_barcodes_reads"] += 1

        aligned_sample = aligned_by_sample.get(sample, {})
        nodup_sample = nodup_by_sample.get(sample, {})
        for group in group_order[sample]:
            for mark in mark_order.get(sample, []):
                key = (group, mark)
                aligned_bam = aligned_sample.get(key)
                nodup_bam = nodup_sample.get(key)
                if aligned_bam is None:
                    raise ValueError(f"Missing DNA aligned BAM for sample '{sample}' group '{group}' mark '{mark}'")
                if nodup_bam is None:
                    raise ValueError(f"Missing DNA NoDup BAM for sample '{sample}' group '{group}' mark '{mark}'")

                aligned_reads = samtools_count(aligned_bam)
                nodup_reads = samtools_count(nodup_bam)

                group_stage_totals[group]["properly_mapped_reads"] += aligned_reads
                group_stage_totals[group]["unique_nodup_reads"] += nodup_reads
                group_mark_stage_totals[key]["properly_mapped_reads"] = aligned_reads
                group_mark_stage_totals[key]["unique_nodup_reads"] = nodup_reads

        for group in group_order[sample]:
            for stage in DNA_STAGE_ORDER:
                group_rows.append(
                    {
                        "sample": sample,
                        "group": group,
                        "stage": stage,
                        "reads": int(group_stage_totals[group][stage]),
                    }
                )

        for group in group_order[sample]:
            for mark in mark_order.get(sample, []):
                for stage in DNA_MARK_STAGE_ORDER:
                    group_mark_rows.append(
                        {
                            "sample": sample,
                            "group": group,
                            "mark": mark,
                            "stage": stage,
                            "reads": int(group_mark_stage_totals[(group, mark)][stage]),
                        }
                    )

        for stage in DNA_STAGE_ORDER:
            sample_rows.append(
                {
                    "sample": sample,
                    "stage": stage,
                    "reads": int(sum(group_stage_totals[group][stage] for group in group_order[sample])),
                }
            )

    group_df = pd.DataFrame(group_rows)
    group_mark_df = pd.DataFrame(group_mark_rows)
    sample_df = pd.DataFrame(sample_rows)
    write_table(group_df, outdir / "dna_group_stage_counts.tsv")
    write_table(group_mark_df, outdir / "dna_group_mark_stage_counts.tsv")
    write_table(sample_df, outdir / "dna_sample_stage_counts.tsv")

    for sample in sample_ids:
        plot_group_overview(
            group_df,
            sample,
            "group",
            DNA_STAGE_ORDER,
            dna_plot_dir / f"{sample}.group_stage_counts.png",
            f"DNA QC: {sample}",
        )
        for group in group_order[sample]:
            plot_single_group(
                group_df,
                sample,
                group,
                DNA_STAGE_ORDER,
                group_plot_dir / f"{sample}__{group}.png",
                f"DNA QC: {sample} / {group}",
            )
            plot_group_marks(
                group_mark_df,
                sample,
                group,
                DNA_MARK_STAGE_ORDER,
                group_mark_plot_dir / f"{sample}__{group}.png",
                f"DNA QC: {sample} / {group} by mark",
            )

    return 0


def main() -> int:
    args = parse_args()
    if args.mode == "rna":
        return build_rna_qc(args)
    if args.mode == "dna":
        return build_dna_qc(args)
    raise ValueError(f"Unsupported mode: {args.mode}")


if __name__ == "__main__":
    sys.exit(main())
