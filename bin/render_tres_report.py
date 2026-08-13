#!/usr/bin/env python3
"""Render a compact TrESFlow HTML QC report from pipeline artifacts."""

from __future__ import annotations

import argparse
import html
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean


def parse_number(value: str):
    value = value.strip().replace(",", "")
    if value in {"", "-nan", "nan", "NaN", "N/A"}:
        return None
    try:
        if re.match(r"^-?\d+$", value):
            return int(value)
        return float(value)
    except ValueError:
        return None


def parse_percent(value: str):
    value = value.strip()
    if not value:
        return None
    if value.endswith("%"):
        value = value[:-1]
    parsed = parse_number(value)
    return float(parsed) if parsed is not None else None


def pct_from_fraction(value):
    if value is None:
        return None
    return float(value) * 100.0


def fmt_int(value):
    if value is None:
        return "n/a"
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return "n/a"


def fmt_pct(value):
    if value is None:
        return "n/a"
    return f"{float(value):.1f}%"


def pct(numerator, denominator):
    if numerator is None or denominator in {None, 0}:
        return None
    return float(numerator) / float(denominator) * 100.0


def fmt_count_pct(count, denominator):
    percent = pct(count, denominator)
    if percent is None:
        return "n/a"
    return f"{fmt_pct(percent)} ({fmt_int(count)})"


def css_width(value):
    if value is None:
        return "0"
    return str(max(0.0, min(100.0, float(value))))


def read_tsv_stats(path: Path):
    stats = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        key = parts[0]
        count = parse_number(parts[1]) if len(parts) > 1 else None
        percent = parse_percent(parts[2]) if len(parts) > 2 else None
        stats[key] = {"count": count, "percent": percent, "raw": parts}
    return stats


def read_counts(path: Path):
    total = 0
    observed = 0
    top = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if not parts:
            continue
        count = parse_number(parts[0])
        if count is not None:
            total += int(count)
            observed += 1
            if len(top) < 5:
                top.append({"count": int(count), "value": parts[1] if len(parts) > 1 else ""})
    return {"total": total, "observed": observed, "top": top}


def read_csv_key_values(path: Path):
    values = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or "," not in line:
            continue
        key, value = line.split(",", 1)
        values[key.strip()] = parse_number(value)
    return values


def read_flagstat(path: Path):
    metrics = {"file": path.name}
    for line in path.read_text(errors="replace").splitlines():
        total = re.match(r"^(\d+) \+ \d+ in total", line)
        primary = re.match(r"^(\d+) \+ \d+ primary$", line)
        mapped = re.match(r"^(\d+) \+ \d+ mapped \(([^%]+)%", line)
        primary_mapped = re.match(r"^(\d+) \+ \d+ primary mapped \(([^%]+)%", line)
        paired = re.match(r"^(\d+) \+ \d+ paired in sequencing", line)
        read1 = re.match(r"^(\d+) \+ \d+ read1", line)
        read2 = re.match(r"^(\d+) \+ \d+ read2", line)
        properly = re.match(r"^(\d+) \+ \d+ properly paired \(([^%]+)%", line)
        duplicates = re.match(r"^(\d+) \+ \d+ duplicates", line)
        if total:
            metrics["total"] = int(total.group(1))
        elif primary:
            metrics["primary"] = int(primary.group(1))
        elif mapped:
            metrics["mapped"] = int(mapped.group(1))
            metrics["mapped_percent"] = float(mapped.group(2))
        elif primary_mapped:
            metrics["primary_mapped"] = int(primary_mapped.group(1))
            metrics["primary_mapped_percent"] = float(primary_mapped.group(2))
        elif paired:
            metrics["paired_in_sequencing"] = int(paired.group(1))
        elif read1:
            metrics["read1"] = int(read1.group(1))
        elif read2:
            metrics["read2"] = int(read2.group(1))
        elif properly:
            metrics["properly_paired"] = int(properly.group(1))
            metrics["properly_paired_percent"] = float(properly.group(2))
        elif duplicates:
            metrics["duplicates"] = int(duplicates.group(1))
    return metrics


def read_duplicate_metrics(path: Path):
    header = None
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            if "PERCENT_DUPLICATION" not in parts:
                continue
            header = parts
            continue
        values = dict(zip(header, parts))
        percent_dup = parse_number(values.get("PERCENT_DUPLICATION", ""))
        read_pairs = parse_number(values.get("READ_PAIRS_EXAMINED", ""))
        read_pair_duplicates = parse_number(values.get("READ_PAIR_DUPLICATES", ""))
        optical_duplicates = parse_number(values.get("READ_PAIR_OPTICAL_DUPLICATES", ""))
        estimated_library_size = parse_number(values.get("ESTIMATED_LIBRARY_SIZE", ""))
        unpaired_reads = parse_number(values.get("UNPAIRED_READS_EXAMINED", ""))
        unpaired_duplicates = parse_number(values.get("UNPAIRED_READ_DUPLICATES", ""))
        return {
            "file": path.name,
            "library": values.get("LIBRARY"),
            "percent_duplication": float(percent_dup) * 100.0 if percent_dup is not None else None,
            "unique_percent": (1.0 - float(percent_dup)) * 100.0 if percent_dup is not None else None,
            "read_pairs_examined": int(read_pairs) if read_pairs is not None else None,
            "read_pair_duplicates": int(read_pair_duplicates) if read_pair_duplicates is not None else None,
            "read_pair_optical_duplicates": int(optical_duplicates) if optical_duplicates is not None else None,
            "estimated_library_size": int(estimated_library_size) if estimated_library_size is not None else None,
            "unpaired_reads_examined": int(unpaired_reads) if unpaired_reads is not None else None,
            "unpaired_read_duplicates": int(unpaired_duplicates) if unpaired_duplicates is not None else None,
        }
    return {"file": path.name}


def strip_suffix(name: str, suffix: str) -> str:
    return name[: -len(suffix)] if name.endswith(suffix) else name


def collect_inputs(input_dir: Path):
    files = [p for p in input_dir.rglob("*") if p.is_file()]
    samples = defaultdict(lambda: {"rna": {}, "dna": {}})
    rna_summaries = []
    star_logs = []
    flagstats = []
    duplicate_metrics = []
    samtools_stats = []

    for path in files:
        name = path.name

        if name.endswith(".rna_sample_barcode.stats.tsv"):
            sample = strip_suffix(name, ".rna_sample_barcode.stats.tsv")
            samples[sample]["sample_id"] = sample
            samples[sample]["rna"]["sample_barcode"] = read_tsv_stats(path)
        elif name.endswith(".sample_barcode.stats.tsv") and not name.endswith(".dna_sample_barcode.stats.tsv"):
            sample = strip_suffix(name, ".sample_barcode.stats.tsv")
            samples[sample]["sample_id"] = sample
            samples[sample]["rna"]["sample_barcode"] = read_tsv_stats(path)
        elif name.endswith(".dna_sample_barcode.stats.tsv"):
            sample = strip_suffix(name, ".dna_sample_barcode.stats.tsv")
            samples[sample]["sample_id"] = sample
            samples[sample]["dna"]["sample_barcode"] = read_tsv_stats(path)
        elif name.endswith(".dna_modality.stats.tsv"):
            sample = strip_suffix(name, ".dna_modality.stats.tsv")
            samples[sample]["sample_id"] = sample
            samples[sample]["dna"]["modality_barcode"] = read_tsv_stats(path)
        elif name.endswith(".rna_umi.counts.tsv"):
            sample = strip_suffix(name, ".rna_umi.counts.tsv")
            samples[sample]["sample_id"] = sample
            samples[sample]["rna"]["umi_counts"] = read_counts(path)
        elif name.endswith(".umi.counts.tsv"):
            sample = strip_suffix(name, ".umi.counts.tsv")
            samples[sample]["sample_id"] = sample
            samples[sample]["rna"]["umi_counts"] = read_counts(path)
        elif name.endswith(".rna_cell.stats_L1.tsv") or name.endswith(".rna_cell.stats_L2.tsv") or name.endswith(".rna_cell.stats_L3.tsv"):
            sample = re.sub(r"\.rna_cell\.stats_L[123]\.tsv$", "", name)
            samples[sample]["sample_id"] = sample
            samples[sample]["rna"].setdefault("cell_barcode", []).append(read_tsv_stats(path))
        elif name.endswith(".cell.stats_L1.tsv") or name.endswith(".cell.stats_L2.tsv") or name.endswith(".cell.stats_L3.tsv"):
            sample = re.sub(r"\.cell\.stats_L[123]\.tsv$", "", name)
            samples[sample]["sample_id"] = sample
            samples[sample]["rna"].setdefault("cell_barcode", []).append(read_tsv_stats(path))
        elif name.endswith(".dna_cell.stats_L1.tsv") or name.endswith(".dna_cell.stats_L2.tsv") or name.endswith(".dna_cell.stats_L3.tsv"):
            sample = re.sub(r"\.dna_cell\.stats_L[123]\.tsv$", "", name)
            samples[sample]["sample_id"] = sample
            samples[sample]["dna"].setdefault("cell_barcode", []).append(read_tsv_stats(path))
        elif name == "Summary.csv":
            split = path.parent.name.replace(".Solo.outGeneFull", "")
            values = read_csv_key_values(path)
            values["split_id"] = split
            rna_summaries.append(values)
        elif name.endswith(".Log.final.out"):
            star_logs.append({"file": name})
        elif name.endswith(".flagstat"):
            metrics = read_flagstat(path)
            metrics["id"] = strip_suffix(name, ".flagstat")
            flagstats.append(metrics)
        elif name.endswith(".DuplicateMetrics.txt"):
            duplicate_metrics.append(read_duplicate_metrics(path))
        elif name.endswith(".stats"):
            samtools_stats.append({"file": name})

    return {
        "samples": dict(samples),
        "rna_summaries": rna_summaries,
        "star_logs": star_logs,
        "flagstats": flagstats,
        "duplicate_metrics": duplicate_metrics,
        "samtools_stats": samtools_stats,
        "input_file_count": len(files),
    }


def weighted_rna_percent(records, key):
    numerator = 0.0
    denominator = 0.0
    for record in records:
        reads = record.get("Number of Reads")
        value = record.get(key)
        if reads is None or value is None:
            continue
        numerator += float(reads) * float(value)
        denominator += float(reads)
    if denominator == 0:
        return None
    return pct_from_fraction(numerator / denominator)


def average(values):
    clean = [float(v) for v in values if v is not None]
    return mean(clean) if clean else None


def weighted_dna_unique_percent(records):
    duplicate_total = 0
    examined_total = 0
    for record in records:
        read_pairs = record.get("read_pairs_examined")
        read_pair_duplicates = record.get("read_pair_duplicates")
        unpaired_reads = record.get("unpaired_reads_examined") or 0
        unpaired_duplicates = record.get("unpaired_read_duplicates") or 0
        if read_pairs is None or read_pair_duplicates is None:
            continue
        duplicate_total += int(read_pair_duplicates) + int(unpaired_duplicates)
        examined_total += int(read_pairs) + int(unpaired_reads)
    if examined_total:
        return (1.0 - (duplicate_total / examined_total)) * 100.0
    return average([item.get("unique_percent") for item in records])


def cell_barcode_percent(cell_stats):
    if not cell_stats:
        return None
    candidates = []
    for stats in cell_stats:
        if "reads_with_all_ligations" in stats:
            candidates.append(stats["reads_with_all_ligations"]["percent"])
        elif "reads_with_ligation_all_segments" in stats:
            candidates.append(stats["reads_with_ligation_all_segments"]["percent"])
        elif "reads_with_ligation" in stats:
            candidates.append(stats["reads_with_ligation"]["percent"])
        elif "reads_with_L1" in stats:
            candidates.append(stats["reads_with_L1"]["percent"])
    return average(candidates)


def sample_for_split(split_id: str, sample_ids):
    matches = [sample_id for sample_id in sample_ids if split_id == sample_id or split_id.startswith(f"{sample_id}_")]
    return max(matches, key=len) if matches else split_id.split("_", 1)[0]


def flagstat_read_units(record):
    if record.get("read1") is not None:
        return record.get("read1")
    if record.get("paired_in_sequencing") is not None:
        return int(record["paired_in_sequencing"] / 2)
    return record.get("primary_mapped") or record.get("mapped") or record.get("total")


def split_from_flagstat_id(flagstat_id: str):
    split_id = re.sub(r"^(rna|dna)\.", "", flagstat_id)
    split_id = re.sub(r"\.(filtered_cells|aligned|markeddup|nodup)$", "", split_id)
    return split_id


def add_sample_value(sample_metrics, sample_id, key, value):
    if value is None:
        return
    sample_metrics[sample_id][key] = sample_metrics[sample_id].get(key, 0) + int(round(float(value)))


def rna_mapping_counts(records, sample_ids, key):
    counts = defaultdict(int)
    for record in records:
        split_id = record.get("split_id")
        reads = record.get("Number of Reads")
        fraction = record.get(key)
        if not split_id or reads is None or fraction is None:
            continue
        sample_id = sample_for_split(split_id, sample_ids)
        counts[sample_id] += int(round(float(reads) * float(fraction)))
    return counts


def flagstat_counts_by_sample(flagstats, sample_ids, modality, stage):
    counts = defaultdict(dict)
    needle = f".{stage}"
    for record in flagstats:
        flagstat_id = record.get("id", "")
        if not flagstat_id.startswith(f"{modality}.") or needle not in flagstat_id:
            continue
        sample_id = sample_for_split(split_from_flagstat_id(flagstat_id), sample_ids)
        add_sample_value(counts, sample_id, "reads", flagstat_read_units(record))
    return {sample_id: values.get("reads") for sample_id, values in counts.items()}


def build_metrics(collected, library_name="unknown library"):
    samples = collected["samples"]
    library_name = library_name or "unknown library"
    sample_ids = sorted(samples)
    rna_transcriptome_counts = rna_mapping_counts(
        collected["rna_summaries"],
        sample_ids,
        "Reads Mapped to GeneFull: Unique GeneFull",
    )
    rna_genome_counts = rna_mapping_counts(
        collected["rna_summaries"],
        sample_ids,
        "Reads Mapped to Genome: Unique",
    )
    rna_usable_counts = flagstat_counts_by_sample(collected["flagstats"], sample_ids, "rna", "filtered_cells")
    dna_mapped_counts = flagstat_counts_by_sample(collected["flagstats"], sample_ids, "dna", "aligned")
    dna_usable_counts = flagstat_counts_by_sample(collected["flagstats"], sample_ids, "dna", "nodup")

    sequencing_rows = []
    for sample_id in sorted(samples):
        sample = samples[sample_id]

        rna = sample.get("rna", {})
        rna_sb = rna.get("sample_barcode", {})
        if rna_sb:
            sequencing_rows.append(
                {
                    "library_name": library_name,
                    "sample_id": sample_id,
                    "modality": "RNA",
                    "reads": rna_sb.get("reads", {}).get("count"),
                    "confidently_mapped": pct(rna_genome_counts.get(sample_id), rna_sb.get("reads", {}).get("count")),
                    "confidently_mapped_reads": rna_genome_counts.get(sample_id),
                    "valid_sample_barcodes": rna_sb.get("bc_reads", {}).get("percent"),
                    "valid_cell_barcodes": cell_barcode_percent(rna.get("cell_barcode", [])),
                    "valid_modality_barcodes": None,
                    "usable_reads": rna_usable_counts.get(sample_id),
                }
            )

        dna = sample.get("dna", {})
        dna_sb = dna.get("sample_barcode", {})
        if dna_sb:
            sequencing_rows.append(
                {
                    "library_name": library_name,
                    "sample_id": sample_id,
                    "modality": "DNA",
                    "reads": dna_sb.get("reads", {}).get("count"),
                    "confidently_mapped": pct(dna_mapped_counts.get(sample_id), dna_sb.get("reads", {}).get("count")),
                    "confidently_mapped_reads": dna_mapped_counts.get(sample_id),
                    "valid_sample_barcodes": dna_sb.get("bc_reads", {}).get("percent"),
                    "valid_modality_barcodes": dna.get("modality_barcode", {}).get("bc_reads", {}).get("percent"),
                    "valid_cell_barcodes": cell_barcode_percent(dna.get("cell_barcode", [])),
                    "usable_reads": dna_usable_counts.get(sample_id),
                }
            )

    rna_raw_reads = sum(row.get("reads") or 0 for row in sequencing_rows if row["modality"] == "RNA")
    dna_raw_reads = sum(row.get("reads") or 0 for row in sequencing_rows if row["modality"] == "DNA")
    rna_transcriptome_count = sum(rna_transcriptome_counts.values())
    rna_genome_count = sum(rna_genome_counts.values())
    dna_mapped_count = sum(dna_mapped_counts.values())
    dna_usable_count = sum(dna_usable_counts.values())

    mapping_quality = {
        "rna_confidently_mapped_to_transcriptome_reads": rna_transcriptome_count,
        "rna_confidently_mapped_to_transcriptome_percent": pct(rna_transcriptome_count, rna_raw_reads),
        "rna_confidently_mapped_to_genome_reads": rna_genome_count,
        "rna_confidently_mapped_to_genome_percent": pct(rna_genome_count, rna_raw_reads),
        "dna_confidently_mapped_reads": dna_mapped_count,
        "dna_confidently_mapped_percent": pct(dna_mapped_count, dna_raw_reads),
        "dna_unique_reads": dna_usable_count,
        "dna_unique_reads_percent": pct(dna_usable_count, dna_raw_reads),
    }
    main_statistics = [
        {
            "metric": "RNA confidently mapped to transcriptome",
            "value": mapping_quality["rna_confidently_mapped_to_transcriptome_percent"],
            "count": rna_transcriptome_count,
            "denominator": rna_raw_reads,
            "value_type": "percent",
            "subtitle": "STARsolo GeneFull unique reads / raw RNA reads",
            "modality": "RNA",
        },
        {
            "metric": "DNA unique reads",
            "value": mapping_quality["dna_unique_reads_percent"],
            "count": dna_usable_count,
            "denominator": dna_raw_reads,
            "value_type": "percent",
            "subtitle": "Final NoDup BAM read pairs / raw DNA reads",
            "modality": "DNA",
        },
    ]
    libraries = [
        {
            "library_name": library_name,
            "mapping_quality": mapping_quality,
            "main_statistics": main_statistics,
            "sequencing_quality": sequencing_rows,
        }
    ]
    export_rows = build_export_rows(libraries)

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "library_name": library_name,
        "libraries": libraries,
        "mapping_quality": mapping_quality,
        "sequencing_quality": sequencing_rows,
        "export_rows": export_rows,
        "inputs": {
            "file_count": collected["input_file_count"],
            "rna_summary_count": len(collected["rna_summaries"]),
            "flagstat_count": len(collected["flagstats"]),
            "duplicate_metrics_count": len(collected["duplicate_metrics"]),
            "samtools_stats_count": len(collected["samtools_stats"]),
        },
        "raw": collected,
    }


def export_value(value, value_type=None):
    if value_type == "percent":
        return fmt_pct(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.4g}"
    if value is None:
        return "n/a"
    return str(value)


def build_export_rows(libraries):
    rows = []
    for library in libraries:
        library_name = library["library_name"]
        for card in library["main_statistics"]:
            rows.append(
                {
                    "section": "Main statistics",
                    "library": library_name,
                    "modality": card["modality"],
                    "fastq_id": "",
                    "metric": card["metric"],
                    "value": export_value(card.get("value"), card.get("value_type")),
                    "absolute_reads": fmt_int(card.get("count")),
                    "raw_reads": fmt_int(card.get("denominator")),
                    "number_of_reads": "",
                    "confidently_mapped": "",
                    "valid_sample_barcodes": "",
                    "valid_cell_barcodes": "",
                    "valid_modality_barcodes": "",
                    "usable_reads": "",
                    "details": card["subtitle"],
                }
            )
        for row in library["sequencing_quality"]:
            rows.append(
                {
                    "section": "Sequencing quality",
                    "library": library_name,
                    "modality": row["modality"],
                    "fastq_id": row["sample_id"],
                    "metric": "sample-level sequencing QC",
                    "value": "",
                    "number_of_reads": fmt_int(row.get("reads")),
                    "confidently_mapped": fmt_count_pct(row.get("confidently_mapped_reads"), row.get("reads")),
                    "valid_sample_barcodes": fmt_pct(row.get("valid_sample_barcodes")),
                    "valid_modality_barcodes": fmt_pct(row.get("valid_modality_barcodes"))
                    if row["modality"] == "DNA"
                    else "n/a",
                    "valid_cell_barcodes": fmt_pct(row.get("valid_cell_barcodes")),
                    "usable_reads": fmt_int(row.get("usable_reads")),
                    "details": "",
                }
            )
    return rows


def metric_card(title, value, subtitle, modality=None):
    modality_class = (modality or "").lower()
    absolute_line = ""
    if isinstance(value, dict):
        count = value.get("count")
        denominator = value.get("denominator")
        value = value.get("percent")
        absolute_line = f"<div class=\"metric-count\">{fmt_int(count)} / {fmt_int(denominator)} raw reads</div>"
    return f"""
    <div class="metric-card {html.escape(modality_class)}">
      <div class="metric-title">{html.escape(title)}</div>
      <div class="metric-value">{fmt_pct(value)}</div>
      {absolute_line}
      <div class="bar"><span style="width: {css_width(value)}%"></span></div>
      <div class="metric-subtitle">{html.escape(subtitle)}</div>
    </div>
    """


def sequencing_table(rows):
    row_html = []
    for row in rows:
        modality_class = "rna" if row["modality"] == "RNA" else "dna"
        row_html.append(
            f"""
            <tr>
              <td><span class="pill {modality_class}">{html.escape(row['modality'])}</span></td>
              <td>{html.escape(row['sample_id'])}</td>
              <td class="num">{fmt_int(row.get('reads'))}</td>
              <td class="num">{fmt_pct(row.get('valid_sample_barcodes'))}</td>
              <td class="num">{fmt_pct(row.get('valid_modality_barcodes')) if row['modality'] == 'DNA' else 'n/a'}</td>
              <td class="num">{fmt_pct(row.get('valid_cell_barcodes'))}</td>
              <td class="num">{fmt_count_pct(row.get('confidently_mapped_reads'), row.get('reads'))}</td>
              <td class="num">{fmt_int(row.get('usable_reads'))}</td>
            </tr>
            """
        )

    return f"""
      <table class="detail-table">
        <thead>
          <tr>
            <th>Modality</th>
            <th>Fastq ID</th>
            <th class="num">Number of reads</th>
            <th class="num">Valid sample barcodes</th>
            <th class="num">Valid modality barcodes</th>
            <th class="num">Valid cell barcodes</th>
            <th class="num">Confidently mapped to genome</th>
            <th class="num">Usable reads</th>
          </tr>
        </thead>
        <tbody>
          {''.join(row_html) if row_html else '<tr><td colspan="8">No sequencing-quality inputs were detected.</td></tr>'}
        </tbody>
      </table>
    """


def render_library(library):
    cards = []
    for card in library["main_statistics"]:
        cards.append(
            metric_card(
                card["metric"],
                {
                    "percent": card.get("value"),
                    "count": card.get("count"),
                    "denominator": card.get("denominator"),
                },
                card["subtitle"],
                card.get("modality"),
            )
        )

    return f"""
    <section class="library-block">
      <section class="panel">
        <div class="panel-header">
          <h2>Mapping Quality</h2>
        </div>
        <div class="grid">
          {''.join(cards)}
        </div>
      </section>

      <section class="panel">
        <div class="panel-header">
          <h2>Sequencing Quality</h2>
        </div>
        {sequencing_table(library['sequencing_quality'])}
      </section>
    </section>
    """


def render_html(metrics):
    libraries = metrics.get("libraries") or []
    library_label = metrics.get("library_name") or "unknown library"
    export_rows_json = json.dumps(metrics.get("export_rows", []), ensure_ascii=False).replace("</", "<\\/")

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>TrESFlow QC Report</title>
  <style>
    :root {{
      --ink: #10284a;
      --muted: #5f6f8b;
      --line: #dbe2ea;
      --panel: #ffffff;
      --bg: #f4f7fb;
      --rna: #1b7f6a;
      --dna: #ff3b30;
      --good: #246bfe;
      --accent: #3d1c96;
    }}
    body {{
      margin: 0;
      padding: 28px;
      background: var(--bg);
      color: var(--ink);
      font: 16px/1.5 "Aptos", "Segoe UI", sans-serif;
    }}
    .page {{
      max-width: 1240px;
      margin: 0 auto;
    }}
    .hero {{
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 24px;
      margin-bottom: 24px;
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 34px;
      letter-spacing: -0.03em;
    }}
    .library-title {{
      color: var(--accent);
      margin-left: 12px;
      white-space: nowrap;
    }}
    .meta {{
      color: var(--muted);
      font-size: 14px;
    }}
    .actions {{
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      justify-content: flex-end;
    }}
    button {{
      border: 1px solid #c9d4e3;
      border-radius: 999px;
      background: #fff;
      color: var(--ink);
      cursor: pointer;
      font: inherit;
      font-weight: 700;
      padding: 9px 14px;
      box-shadow: 0 1px 4px rgba(16, 40, 74, 0.08);
    }}
    button.primary {{
      background: var(--ink);
      border-color: var(--ink);
      color: #fff;
    }}
    .library-block {{
      margin: 24px 0 34px;
    }}
    .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 12px;
      box-shadow: 0 2px 8px rgba(16, 40, 74, 0.08);
      margin: 18px 0;
      overflow: hidden;
    }}
    .panel-header {{
      display: flex;
      align-items: baseline;
      padding: 20px 22px;
      border-bottom: 1px solid var(--line);
    }}
    h2 {{
      margin: 0;
      font-size: 25px;
      letter-spacing: -0.02em;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 1px;
      background: var(--line);
    }}
    .metric-card {{
      background: #fff;
      padding: 22px;
      min-height: 150px;
      position: relative;
    }}
    .metric-card.rna::before,
    .metric-card.dna::before {{
      content: "";
      position: absolute;
      top: 0;
      left: 0;
      right: 0;
      height: 4px;
      background: var(--rna);
    }}
    .metric-card.dna::before {{ background: var(--dna); }}
    .metric-title {{
      color: var(--muted);
      font-weight: 700;
      min-height: 48px;
    }}
    .metric-value {{
      margin-top: 14px;
      font-size: 32px;
      font-weight: 800;
      letter-spacing: -0.03em;
    }}
    .metric-count {{
      margin-top: 2px;
      color: var(--ink);
      font-size: 15px;
      font-weight: 700;
    }}
    .metric-subtitle {{
      margin-top: 10px;
      color: var(--muted);
      font-size: 13px;
    }}
    .bar {{
      height: 10px;
      background: #e9eef5;
      border-radius: 99px;
      overflow: hidden;
      margin-top: 12px;
    }}
    .bar span {{
      display: block;
      height: 100%;
      background: linear-gradient(90deg, var(--good), #62d2a2);
      border-radius: 99px;
    }}
    .detail-table {{
      table-layout: fixed;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
    }}
    th, td {{
      border-bottom: 1px solid var(--line);
      padding: 16px 18px;
      vertical-align: middle;
    }}
    th {{
      color: var(--muted);
      text-align: left;
      font-size: 14px;
      letter-spacing: 0.02em;
    }}
    td.num, th.num {{
      text-align: right;
      font-variant-numeric: tabular-nums;
    }}
    .pill {{
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 4px 10px;
      color: white;
      font-size: 12px;
      font-weight: 800;
      letter-spacing: 0.04em;
    }}
    .pill.rna {{ background: var(--rna); }}
    .pill.dna {{ background: var(--dna); }}
    .empty {{
      color: var(--muted);
      padding: 22px;
    }}
    @media (max-width: 900px) {{
      body {{ padding: 14px; }}
      .grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      .panel-header {{ display: block; }}
      .hero {{ display: block; }}
      .actions {{ justify-content: flex-start; margin-top: 14px; }}
      .library-title {{ display: block; margin-left: 0; }}
    }}
    @media (max-width: 620px) {{
      .grid {{ grid-template-columns: 1fr; }}
      table {{ font-size: 14px; }}
      th, td {{ padding: 10px; }}
    }}
  </style>
</head>
<body>
  <main class="page">
    <section class="hero">
      <div>
        <h1>TrESFlow QC Report <span class="library-title">{html.escape(library_label)}</span></h1>
        <div class="meta">Generated {html.escape(metrics['generated_at'])}</div>
      </div>
      <div class="actions">
        <button class="primary" type="button" onclick="downloadCsv()">Export CSV</button>
        <button type="button" onclick="downloadExcel()">Export Excel</button>
      </div>
    </section>

    {''.join(render_library(library) for library in libraries) if libraries else '<div class="panel empty">No library metrics were detected.</div>'}
  </main>
  <script id="export-data" type="application/json">{export_rows_json}</script>
  <script>
    const exportRows = JSON.parse(document.getElementById("export-data").textContent || "[]");

    function exportHeaders(rows) {{
      const preferred = [
        "section",
        "library",
        "modality",
        "fastq_id",
        "metric",
        "value",
        "absolute_reads",
        "raw_reads",
        "number_of_reads",
        "confidently_mapped",
        "valid_sample_barcodes",
        "valid_modality_barcodes",
        "valid_cell_barcodes",
        "usable_reads",
        "details"
      ];
      const extra = Array.from(new Set(rows.flatMap(row => Object.keys(row)))).filter(key => !preferred.includes(key));
      return preferred.concat(extra);
    }}

    function csvEscape(value) {{
      const text = value === null || value === undefined ? "" : String(value);
      return /[",\\n\\r]/.test(text) ? '"' + text.replace(/"/g, '""') + '"' : text;
    }}

    function downloadBlob(filename, content, type) {{
      const options = Object.create(null);
      options.type = type;
      const blob = new Blob([content], options);
      const link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      URL.revokeObjectURL(link.href);
      link.remove();
    }}

    function downloadCsv() {{
      const headers = exportHeaders(exportRows);
      const lines = [headers.join(",")].concat(
        exportRows.map(row => headers.map(header => csvEscape(row[header])).join(","))
      );
      downloadBlob("tresflow_qc_report.csv", lines.join("\\n") + "\\n", "text/csv;charset=utf-8");
    }}

    function htmlEscape(value) {{
      return String(value === null || value === undefined ? "" : value)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;");
    }}

    function downloadExcel() {{
      const headers = exportHeaders(exportRows);
      const headerHtml = headers.map(header => "<th>" + htmlEscape(header) + "</th>").join("");
      const bodyHtml = exportRows.map(row =>
        "<tr>" + headers.map(header => "<td>" + htmlEscape(row[header]) + "</td>").join("") + "</tr>"
      ).join("");
      const workbook = '<html><head><meta charset="utf-8"></head><body><table>' +
        "<thead><tr>" + headerHtml + "</tr></thead><tbody>" + bodyHtml + "</tbody></table></body></html>";
      downloadBlob("tresflow_qc_report.xls", workbook, "application/vnd.ms-excel;charset=utf-8");
    }}
  </script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-html", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--library-name", default="unknown library")
    args = parser.parse_args()

    collected = collect_inputs(args.input_dir)
    metrics = build_metrics(collected, args.library_name)

    args.output_json.write_text(json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8")
    args.output_html.write_text(render_html(metrics), encoding="utf-8")


if __name__ == "__main__":
    main()
