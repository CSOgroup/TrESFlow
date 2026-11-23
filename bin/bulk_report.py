#!/usr/bin/env python3
"""Render pilot QC from declared process outputs and combine raw integer RNA counts."""
import base64
import csv
import json
import sys
from html import escape
from itertools import groupby, zip_longest
from contextlib import ExitStack
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from tresflow_qc.bulk_barcodes import abundance_svg, label, raw_panels, routed_panels
from tresflow_qc.core import read_barcode_counts
from tresflow_qc.html_report import CSS, PLOT_EXPORT_SCRIPT, plot_markup, prepare_plots
from tresflow_qc.plots import DNA, RNA, _ellipsize_pixels, svg_figure


def merge_counts(paths, output):
    # Every group uses the same annotated STAR reference. Preserve its row order
    # and require an exact feature match before publishing a multi-group matrix.
    with ExitStack() as stack, output.open("w") as out:
        readers = [csv.reader(stack.enter_context(p.open()), delimiter="\t") for p in sorted(paths)]
        for rows in zip_longest(*readers):
            if any(row is None or len(row) != 2 or row[0] != rows[0][0] for row in rows):
                raise ValueError("RNA count files disagree on reference gene identifiers or order")
            print(rows[0][0], *(row[1] for row in rows), sep="\t", file=out)


def number(value):
    if value is None:
        return "unavailable"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def table(headers, rows):
    head = "".join(f"<th>{escape(col)}</th>" for col in headers)
    body = "".join("<tr>" + "".join(f"<td>{escape(number(v))}</td>" for v in row) + "</tr>" for row in rows)
    return f'<div class="table-scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def metric(label, value, denominator=None, scope=""):
    pct = f"{100 * value / denominator:.2f}%" if denominator and value is not None else ("unavailable" if scope else "—")
    return label, value, pct, scope


def yield_svg(points, title, scope, xlabel, ylabel, color, bin_size=0):
    left, right, top, bottom = 65, 430, 85, 178
    xmax = max([x + bin_size for x, _ in points] + [1])
    ymax = max([y for _, y in points] + [1])
    coords = [(left + x*(right-left)/xmax, bottom - y*(bottom-top)/ymax) for x, y in points]
    if bin_size:
        marks = "".join(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bin_size*(right-left)/xmax:.2f}" '
                        f'height="{bottom-y:.2f}" fill="{color}"/>' for x, y in coords)
    else:
        line = " ".join(f"{x:.2f},{y:.2f}" for x, y in coords)
        marks = f'<polyline points="{line}" fill="none" stroke="{color}" stroke-width="2.5"/>'
    body = (f'<text class="panel-title" x="10" y="22">{escape(title)}</text>'
            f'<text class="axis muted" x="10" y="44">{escape(_ellipsize_pixels(scope, 435, 13))}</text>'
            f'<text class="axis" x="{left}" y="70">{escape(ylabel)}</text>'
            f'<path d="M{left} {top}V{bottom}H{right}" fill="none" stroke="#aeb4b4"/>{marks}'
            f'<text class="axis" x="{left-5}" y="{top+5}" text-anchor="end">{ymax:,.0f}</text>'
            f'<text class="axis" x="{left}" y="198">0</text>'
            f'<text class="axis" x="{right}" y="198" text-anchor="end">{xmax:,.0f}</text>'
            f'<text class="axis" x="240" y="217" text-anchor="middle">{escape(xlabel)}</text>')
    return svg_figure(460, 228, body, f"{title} · {scope}")


def rna_metrics(q):
    retained, routed = q["retained_pairs"], q["routed_pairs"]
    mapped = q["genomeU"]
    sense, anti = q.get("sense", 0), q.get("antisense", 0)
    return [
        metric("Routed input pairs", routed),
        metric("Pairs after clipping and length gate", retained, routed, "routed pairs"),
        metric("Uniquely mapped pairs", mapped, retained, "pairs submitted to STAR"),
        metric("Multimapped pairs", q["genomeM"], retained, "pairs submitted to STAR"),
        metric("Uniquely assigned GeneFull pairs", q["featureU"], retained, "pairs submitted to STAR"),
        *(metric(f"{region.capitalize()} pairs", q.get(region, 0), mapped, "uniquely mapped pairs")
          for region in ("exonic", "intronic", "intergenic")),
        metric("Sense pairs", sense, sense + anti, "unambiguous annotated strand"),
        metric("Antisense pairs", anti, sense + anti, "unambiguous annotated strand"),
        metric("Strand ambiguous pairs", q.get("ambiguous", 0)),
        metric("Corrected molecules", q["molecules"], retained, "pairs submitted to STAR"),
        metric("Representatives with one mate mapped", q.get("one_mate_mapped_representatives", 0),
               q["molecules"], "molecule representatives"),
        metric("Pools with molecules", q["pools_detected"]), metric("Detected genes", q["genes_detected"]),
    ]


def dna_metrics(q):
    n, aligned = q["usable_fragments"], q["bwa_primary_pairs"]
    post, mapped, proper = (q[k] for k in ("post_blacklist_primary_pairs", "post_blacklist_mapped_primary_pairs",
                                          "proper_pair_primary_pairs"))
    return [
        metric("BWA input pairs", aligned),
        metric("Mapped primary R1", q.get("bwa_mapped_primary_pairs"), aligned, "input pairs"),
        metric("R1 removed by blacklist", aligned - post, aligned, "input pairs"),
        metric("Unmapped R1 after blacklist", post - mapped, post, "post-blacklist R1"),
        metric("R1 removed by proper-pair filter", mapped - proper, mapped, "post-blacklist mapped R1"),
        metric("Proper-pair R1 outside usable complete fragments", proper - n, proper, "proper-pair R1"),
        metric("Usable fragments", n, aligned, "input pairs"),
        metric("Duplicate-flagged usable fragments", q["usable_duplicate_fragments"], n, "usable fragments"),
        metric("Picard pairs examined", q["read_pairs_examined"]),
        metric("Picard duplicate pairs", q["read_pair_duplicates"], q["read_pairs_examined"], "Picard pairs examined"),
        metric("Picard optical duplicate pairs", q["read_pair_optical_duplicates"], q["read_pair_duplicates"], "duplicate pairs"),
        metric("Peaks", q["peak_count"]), metric("Peak union (bp)", q["peak_bases"]),
        metric("Fragments in peaks (FRiP)", q["fragments_in_peaks"], n, "usable fragments"),
    ] + ([metric(f"Fragments in reference BED ({q['reference_bed']})", q["fragments_in_reference"], n, "usable fragments")]
         if "fragments_in_reference" in q else [])


def records(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


METRIC_HEADERS = ["Metric", "Count", "Percent", "Denominator"]


def details(title, content):
    return f'<details><summary>{escape(title)}</summary>{content}</details>'


def cards(figures):
    return '<div class="bulk-grid">' + "".join(
        f'<article class="plot-card">{plot_markup(item)}</article>' for item in prepare_plots(figures)
    ) + '</div>'


def preprocessing(root):
    panels, tables, warnings = [], [], []
    for path in sorted(root.glob("*_barcode_composition.tsv")):
        rows = records(path)
        sample, mod = rows[0]["sample_id"], rows[0]["modality"]
        prefix = f"{sample}.{mod}"
        scope = f"{sample} · {mod.upper()}"
        with (root / f"{prefix}_sample_barcode.stats.tsv").open() as handle:
            stats = {r[0]: int(r[1]) for r in csv.reader(handle, delimiter="\t")}
        raw, no_sb = stats["reads"], stats["reads_without_bc"]
        gates = [metric("Raw input pairs", raw), metric("Unassigned at sample-barcode tagging", no_sb, raw, "raw input pairs")]
        gates.extend(metric(r["metric"].replace("_", " "), int(r["pairs"]), raw, "raw input pairs")
                     for r in records(root / f"{prefix}_barcode_gates.tsv"))
        tables.append(details(f"{scope} · preprocessing counts", table(METRIC_HEADERS, gates)))
        tag_counts = read_barcode_counts(root / f"{prefix}_sample_barcode.counts.tsv")
        panels.extend(raw_panels(sample, mod, rows, tag_counts, raw, no_sb))
        if no_sb:
            warnings.append(f"{scope} · {no_sb:,} of {raw:,} raw pairs had no accepted sample barcode at tagging.")
        for r in rows:
            if r["barcode_type"] == "sample_barcode" and r["category"] == "NoMatch" and int(r["count"]):
                warnings.append(f"{scope} · {int(r['count']):,} of {int(r['denominator_count']):,} "
                                "ligation-accepted split-input pairs had no sample-barcode group.")
    return panels, tables, warnings


def barcode_section(panels):
    sections = []
    for scope, subset in groupby(sorted(panels, key=lambda p: p.scope), key=lambda p: p.scope):
        figures, exact = [], []
        for p in subset:
            figures.append(dict(kind=f"{p.stage}-{p.tag}", run=scope, svg=abundance_svg(p)))
            for identity, count in sorted(p.counts.items(), key=lambda item: (-item[1], item[0])):
                pct = f"{100 * count / p.denominator:.2f}%" if p.denominator else "unavailable"
                exact.append([p.stage, p.title, label(identity), count, pct, p.denominator, p.definition])
        sections.append(f'<h3>{escape(scope)}</h3>' + cards(figures)
                        + details("Exact barcode counts and denominators", table(
                            ["Stage", "Barcode", "Sequence / group or mark", "Read pairs", "Percent",
                             "Denominator (pairs)", "Denominator definition"], exact)))
    return "".join(sections)


def target_card(title, q, rna):
    if q is None:
        return f'<article class="run-card"><h3>{escape(title)}</h3><p>No routed analysis output was observed.</p></article>'
    if rna:
        fields = [("Routed input pairs", q["routed_pairs"]), ("Corrected molecules", q["molecules"]),
                  ("Genes detected", q["genes_detected"]), ("Pools with molecules", q["pools_detected"])]
    else:
        frip = q["frip"]
        fields = [("BWA input pairs", q["bwa_primary_pairs"]), ("Usable fragments", q["usable_fragments"]),
                  ("Peaks", q["peak_count"]), ("FRiP · usable fragments", "unavailable" if frip is None else f"{100*frip:.2f}%")]
    values = "".join(f"<div><span>{escape(k)}</span><strong>{escape(number(v))}</strong></div>" for k, v in fields)
    return f'<article class="run-card"><h3>{escape(title)}</h3><div class="mini-grid">{values}</div></article>'


def render(spec, root):
    restart = spec["restart"]
    panels, metrics, warnings = ([], [], []) if restart else preprocessing(root)
    targets = sorted(spec["targets"], key=lambda t: (t["sample"], t["group"], t["modality"] != "rna", t["id"]))
    overview, yields, routed = [], [], {}
    for t in targets:
        mod, split = t["modality"], t["id"]
        rna = mod == "rna"
        scope = " · ".join([t["sample"], t["group"], mod.upper()] + ([] if rna else [t["mark"]]))
        path = root / f"{split}.{mod}.qc.json"
        q = json.loads(path.read_text()) if path.exists() else None
        overview.append(target_card(scope, q, rna))
        if q is None:
            continue
        if restart:
            routed[mod, split] = q["barcodes"] if rna else json.loads((root / f"{split}.dna.routed.json").read_text())
        figures = []
        if rna:
            labels = ["Corrected molecules", "Detected genes"]
            note = "Region fractions use aligned blocks of uniquely mapped pairs. Exon overlap takes precedence over gene-body overlap. Strand compares R1 with the unambiguous annotation."
        else:
            bins = {}
            for length, count in q["fragment_lengths"].items():
                lo = int(length) // 25 * 25
                bins[lo] = bins.get(lo, 0) + count
            figures.append(dict(kind="Fragment lengths", run=scope, svg=yield_svg(sorted(bins.items()), "Fragment lengths",
                                scope, "Fragment length (bp, 25 bp bins)", "Usable fragments", DNA, bin_size=25)))
            labels = ["Usable fragments"]
            note = f"Epic2 · {q['caller_status']}. CPM coverage in 50 bp bins includes duplicate-flagged target fragments."
        for i, name in enumerate(labels, 1):
            points = [(0, 0)] + [(row[0], row[i]) for row in q["depth"]]
            figures.append(dict(kind=f"Depth and {name}", run=scope, svg=yield_svg(points, f"{name} by depth", scope,
                                "Retained sequencing pairs" if rna else "Sequencing pairs at BWA",
                                name, RNA if rna else DNA)))
        yields.append(f"<h3>{escape(scope)}</h3>" + cards(figures))
        metrics.append(details(scope + " · analysis metrics", table(METRIC_HEADERS, rna_metrics(q) if rna else dna_metrics(q))
                               + f"<p>{escape(note)}</p>"))

    if restart:
        panels = routed_panels(targets, routed)
    library = ", ".join(sorted({t["library"] for t in targets}))
    stage = ("Measurements start at the supplied routed pairs. Raw depth, upstream losses, unassigned totals "
             "and percentages of raw input are unavailable." if restart else
             "Measurements start at raw input. Sample-barcode tagging and split-input plots have separate denominators.")
    warnings_html = "".join(f'<li>{escape(w)}</li>' for w in warnings)
    if warnings:
        warnings_html = '<aside class="bulk-warnings"><strong>Unassigned sample barcodes</strong><ul>' + warnings_html + '</ul></aside>'
    css = """
    .bulk-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(340px,100%),1fr));gap:.8rem;margin:.75rem 0 1.2rem}
    .bulk-grid .plot-card{min-width:0;padding:.65rem}.bulk-grid .svg-scroll svg{width:100%;max-width:560px;height:auto}
    .bulk-grid .plot-actions{margin-bottom:.2rem}.bulk-warnings{border-left:3px solid var(--orange);padding:.5rem 1rem;background:#fff6ed;margin:1rem 0}
    .bulk-warnings ul{margin:.25rem 0}section{scroll-margin-top:3rem}.bulk-grid+h3{margin-top:1.6rem}
    """
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escape(library)} · Bulk TrES-seq QC</title><style>{CSS}{css}</style></head>
<body><header class="hero"><div class="brand">TrESFlow</div><h1>Bulk TrES-seq QC</h1></header>
<nav aria-label="Report sections"><a href="#overview">Overview</a><a href="#barcodes">Barcode abundance</a>
<a href="#yield">Yield and depth</a><a href="#metrics">Detailed metrics</a></nav><main>
<section id="overview"><h2>{escape(library)}</h2><p class="subtle">{stage}</p>{warnings_html}
<div class="run-cards">{"".join(overview)}</div></section>
<section id="barcodes"><h2>Barcode abundance</h2>{barcode_section(panels)}</section>
<section id="yield"><h2>Yield and sequencing depth</h2><p>Expected yield under independent thinning of observed pairs at 10, 25, 50, 75 and 100% depth. RNA curves use support per corrected molecule and gene with assignments fixed. DNA fragment yield scales linearly because duplicate-flagged fragments are retained.</p>{"".join(yields)}</section>
<section id="metrics"><h2>Detailed metrics</h2>{"".join(metrics)}</section></main>{PLOT_EXPORT_SCRIPT}</body></html>'''

if __name__ == "__main__":
    spec = json.loads(base64.b64decode(sys.argv[1]))
    root = Path(sys.argv[2])
    paths = [root / f"{target['id']}.counts.tsv" for target in spec["targets"] if target["modality"] == "rna"]
    paths = [path for path in paths if path.is_file()]
    if paths:
        merge_counts(paths, Path("gene_by_group.tsv"))
    Path("bulk_report.html").write_text(render(spec, root))
