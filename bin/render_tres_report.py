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
        mapped = re.match(r"^(\d+) \+ \d+ mapped \(([^%]+)%", line)
        properly = re.match(r"^(\d+) \+ \d+ properly paired \(([^%]+)%", line)
        duplicates = re.match(r"^(\d+) \+ \d+ duplicates", line)
        if total:
            metrics["total"] = int(total.group(1))
        elif mapped:
            metrics["mapped"] = int(mapped.group(1))
            metrics["mapped_percent"] = float(mapped.group(2))
        elif properly:
            metrics["properly_paired"] = int(properly.group(1))
            metrics["properly_paired_percent"] = float(properly.group(2))
        elif duplicates:
            metrics["duplicates"] = int(duplicates.group(1))
    return metrics


def read_duplicate_metrics(path: Path):
    header = None
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("##"):
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            continue
        values = dict(zip(header, parts))
        percent_dup = parse_number(values.get("PERCENT_DUPLICATION", ""))
        return {
            "file": path.name,
            "library": values.get("LIBRARY"),
            "percent_duplication": float(percent_dup) * 100.0 if percent_dup is not None else None,
            "unique_percent": (1.0 - float(percent_dup)) * 100.0 if percent_dup is not None else None,
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


def build_metrics(collected):
    samples = collected["samples"]
    rna_transcriptome = weighted_rna_percent(
        collected["rna_summaries"],
        "Reads Mapped to GeneFull: Unique GeneFull",
    )
    rna_genome = weighted_rna_percent(
        collected["rna_summaries"],
        "Reads Mapped to Genome: Unique",
    )

    dna_aligned_flagstats = [
        item for item in collected["flagstats"] if item.get("id", "").startswith("dna.") and ".aligned" in item.get("id", "")
    ]
    dna_mapped = average([item.get("mapped_percent") for item in dna_aligned_flagstats])
    dna_unique = average([item.get("unique_percent") for item in collected["duplicate_metrics"]])

    sequencing_rows = []
    for sample_id in sorted(samples):
        sample = samples[sample_id]

        rna = sample.get("rna", {})
        rna_sb = rna.get("sample_barcode", {})
        if rna_sb:
            sequencing_rows.append(
                {
                    "sample_id": sample_id,
                    "modality": "RNA",
                    "reads": rna_sb.get("reads", {}).get("count"),
                    "valid_sample_barcodes": rna_sb.get("bc_reads", {}).get("percent"),
                    "valid_cell_barcodes": cell_barcode_percent(rna.get("cell_barcode", [])),
                    "valid_umis": None,
                    "observed_umis": rna.get("umi_counts", {}).get("observed"),
                }
            )

        dna = sample.get("dna", {})
        dna_sb = dna.get("sample_barcode", {})
        if dna_sb:
            sequencing_rows.append(
                {
                    "sample_id": sample_id,
                    "modality": "DNA",
                    "reads": dna_sb.get("reads", {}).get("count"),
                    "valid_sample_barcodes": dna_sb.get("bc_reads", {}).get("percent"),
                    "valid_modality_barcodes": dna.get("modality_barcode", {}).get("bc_reads", {}).get("percent"),
                    "valid_cell_barcodes": cell_barcode_percent(dna.get("cell_barcode", [])),
                    "valid_umis": None,
                    "observed_umis": None,
                }
            )

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mapping_quality": {
            "rna_confidently_mapped_to_transcriptome_percent": rna_transcriptome,
            "rna_confidently_mapped_to_genome_percent": rna_genome,
            "dna_confidently_mapped_percent": dna_mapped,
            "dna_unique_reads_percent": dna_unique,
        },
        "sequencing_quality": sequencing_rows,
        "inputs": {
            "file_count": collected["input_file_count"],
            "rna_summary_count": len(collected["rna_summaries"]),
            "flagstat_count": len(collected["flagstats"]),
            "duplicate_metrics_count": len(collected["duplicate_metrics"]),
            "samtools_stats_count": len(collected["samtools_stats"]),
        },
        "raw": collected,
    }


def metric_card(title, value, subtitle):
    return f"""
    <div class="metric-card">
      <div class="metric-title">{html.escape(title)}</div>
      <div class="metric-value">{fmt_pct(value)}</div>
      <div class="bar"><span style="width: {css_width(value)}%"></span></div>
      <div class="metric-subtitle">{html.escape(subtitle)}</div>
    </div>
    """


def render_html(metrics):
    mapping = metrics["mapping_quality"]
    rows = metrics["sequencing_quality"]

    row_html = []
    for row in rows:
        umi_text = "no UMI" if row["modality"] == "DNA" else f"observed {fmt_int(row.get('observed_umis'))}"
        modality_class = "rna" if row["modality"] == "RNA" else "dna"
        row_html.append(
            f"""
            <tr>
              <td><span class="pill {modality_class}">{html.escape(row['modality'])}</span></td>
              <td>{html.escape(row['sample_id'])}</td>
              <td class="num">{fmt_int(row.get('reads'))}</td>
              <td class="num">{fmt_pct(row.get('valid_sample_barcodes'))}</td>
              <td class="num">{fmt_pct(row.get('valid_cell_barcodes'))}</td>
              <td class="num">{fmt_pct(row.get('valid_modality_barcodes')) if row['modality'] == 'DNA' else 'n/a'}</td>
              <td class="num">{html.escape(umi_text)}</td>
            </tr>
            """
        )

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
      margin-bottom: 22px;
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 34px;
      letter-spacing: -0.03em;
    }}
    .meta {{
      color: var(--muted);
      font-size: 14px;
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
      justify-content: space-between;
      gap: 18px;
      padding: 20px 22px;
      border-bottom: 1px solid var(--line);
    }}
    h2 {{
      margin: 0;
      font-size: 25px;
      letter-spacing: -0.02em;
    }}
    .hint {{
      color: var(--muted);
      font-size: 14px;
      text-align: right;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 1px;
      background: var(--line);
    }}
    .metric-card {{
      background: #fff;
      padding: 22px;
      min-height: 150px;
    }}
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
    .notes {{
      padding: 18px 22px;
      color: var(--muted);
      font-size: 14px;
    }}
    @media (max-width: 900px) {{
      body {{ padding: 14px; }}
      .grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      .panel-header {{ display: block; }}
      .hint {{ text-align: left; margin-top: 8px; }}
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
      <h1>TrESFlow QC Report</h1>
      <div class="meta">Generated {html.escape(metrics['generated_at'])}</div>
    </section>

    <section class="panel">
      <div class="panel-header">
        <h2>Mapping Quality</h2>
        <div class="hint">RNA uses STARsolo GeneFull/Genome summary metrics. DNA uses BAM-level samtools/GATK metrics.</div>
      </div>
      <div class="grid">
        {metric_card("RNA confidently mapped to transcriptome", mapping.get("rna_confidently_mapped_to_transcriptome_percent"), "STARsolo GeneFull unique reads")}
        {metric_card("RNA confidently mapped to genome", mapping.get("rna_confidently_mapped_to_genome_percent"), "STARsolo genome unique reads")}
        {metric_card("DNA confidently mapped", mapping.get("dna_confidently_mapped_percent"), "samtools flagstat on aligned DNA BAMs")}
        {metric_card("DNA unique reads", mapping.get("dna_unique_reads_percent"), "1 - GATK MarkDuplicates duplication rate")}
      </div>
    </section>

    <section class="panel">
      <div class="panel-header">
        <h2>Sequencing Quality</h2>
        <div class="hint">DNA has no UMI column by design; RNA UMI is reported as observed extracted UMIs.</div>
      </div>
      <table>
        <thead>
          <tr>
            <th>Modality</th>
            <th>Fastq ID</th>
            <th class="num">Number of reads</th>
            <th class="num">Valid sample barcodes</th>
            <th class="num">Valid cell barcodes</th>
            <th class="num">Valid modality barcodes</th>
            <th class="num">UMI</th>
          </tr>
        </thead>
        <tbody>
          {''.join(row_html) if row_html else '<tr><td colspan="7">No sequencing-quality inputs were detected.</td></tr>'}
        </tbody>
      </table>
      <div class="notes">
        This report is generated from existing TrESFlow artifacts and does not alter pipeline results. The companion JSON file contains the parsed metrics used here.
      </div>
    </section>
  </main>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-html", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()

    collected = collect_inputs(args.input_dir)
    metrics = build_metrics(collected)

    args.output_json.write_text(json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8")
    args.output_html.write_text(render_html(metrics), encoding="utf-8")


if __name__ == "__main__":
    main()
