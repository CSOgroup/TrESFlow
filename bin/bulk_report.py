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
from tresflow_qc.core import read_barcode_counts, read_duplicate_metrics
from tresflow_qc.bulk_fragments import fragment_figure
from tresflow_qc.html_report import CSS, PLOT_EXPORT_SCRIPT, plot_markup, prepare_plots
from tresflow_qc.plots import _complexity_figure, _rna_complexity_figure, mark_palette


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
        metric("Detected genes", q["genes_detected"]),
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


def cards(figures, wide=False):
    return ('<div class="run-plots">' if wide else '<div class="bulk-grid">') + "".join(
        f'<article class="plot-card">{plot_markup(item)}</article>' for item in prepare_plots(figures)
    ) + '</div>'


def percent(value, denominator):
    return f"{100 * value / denominator:.2f}%" if denominator and value is not None else "unavailable"


def saturation(q):
    assigned = q["featureU"]
    return 100 * (1 - q["molecules"] / assigned) if assigned else None


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
            warnings.append(f"{scope} · {percent(no_sb, raw)} of raw input pairs had no accepted sample barcode at tagging.")
        for r in rows:
            if r["barcode_type"] == "sample_barcode" and r["category"] == "NoMatch" and int(r["count"]):
                warnings.append(f"{scope} · {percent(int(r['count']), int(r['denominator_count']))} of "
                                "ligation-accepted pairs had no sample-barcode group.")
    panels = [p for p in panels if p.tag != "SB" or p.stage == "Sample-barcode tagging"]
    return panels, tables, warnings


def barcode_section(panels):
    visible = [p for p in panels if p.tag not in {"L1", "L2", "L3"}]
    ligation = [p for p in panels if p.tag in {"L1", "L2", "L3"}]
    return barcode_panels(visible) + (details("Ligation barcode abundance · L1 / L2 / L3", barcode_panels(ligation)) if ligation else "")


def barcode_panels(panels):
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


def target_card(title, q, rna, raw=None):
    if q is None:
        return f'<article class="run-card"><h3>{escape(title)}</h3><p>No routed analysis output was observed.</p></article>'
    final = q["molecules"] if rna else q["usable_fragments"]
    supplied = q["routed_pairs"] if rna else q["bwa_primary_pairs"]
    denominator = raw if raw is not None else supplied
    source = "raw modality pairs" if raw is not None else "supplied group pairs"
    fields = [
        ("Raw input pairs · modality", number(raw)),
        ("Final molecules" if rna else "Final aligned fragments",
         f"{number(final)} ({percent(final, denominator)} of {source})"),
    ]
    if raw is None:
        fields.append(("Supplied group pairs", number(supplied)))
    if rna:
        value = saturation(q)
        fields.extend([("Genes detected", q["genes_detected"]),
                       ("RNA sequencing saturation", "unavailable" if value is None else f"{value:.2f}%")])
    else:
        fields.extend([
            ("DNA duplicate pairs · Picard", percent(q["read_pair_duplicates"], q["read_pairs_examined"])),
            ("Peaks", q["peak_count"]),
            ("FRiP · final aligned fragments", percent(q["fragments_in_peaks"], final)),
        ])
    values = "".join(f"<div><span>{escape(k)}</span><strong>{escape(number(v))}</strong></div>" for k, v in fields)
    return f'<article class="run-card"><h3>{escape(title)}</h3><div class="mini-grid">{values}</div></article>'


def dna_complexity(q, root, split, scope, branch):
    picard = q.get("picard")
    if picard is None:
        # Older QC JSON can be rendered alongside its published Picard metrics.
        path = root / f"{split}.DuplicateMetrics.txt"
        picard = read_duplicate_metrics(path) if path.exists() else q
    examined = picard.get("read_pairs_examined")
    duplicates = picard.get("read_pair_duplicates")
    optical = picard.get("read_pair_optical_duplicates")
    if any(value is None for value in (examined, duplicates, optical)) or not examined:
        return '<p class="unavailable">Picard duplicate composition and depth prediction unavailable.</p>'
    unique = examined - duplicates
    row = dict(branch=branch, read_pairs_examined=examined, observed_unique_pairs=unique,
               pcr_or_library_duplicate_pairs=duplicates-optical, optical_duplicate_pairs=optical,
               estimated_library_size=picard.get("estimated_library_size"))
    roi = picard.get("roi", [])
    rows = [dict(row, roi_depth_multiplier=p["depth_multiplier"],
                 roi_estimated_unique_pairs=unique*p["yield_multiplier"]) for p in roi] or [row]
    figure = _complexity_figure(scope, rows, mark_palette([branch]))
    unavailable = "" if roi else '<p class="unavailable">Predicted unique pairs with deeper sequencing unavailable · Picard supplied no ROI histogram.</p>'
    return cards([figure], wide=True) + unavailable


def render(spec, root):
    restart = spec["restart"]
    panels, metrics, warnings = ([], [], []) if restart else preprocessing(root)
    targets = sorted(spec["targets"], key=lambda t: (t["sample"], t["group"], t["modality"] != "rna", t["id"]))
    raw_counts = {}
    if not restart:
        for path in root.glob("*_sample_barcode.stats.tsv"):
            sample_mod = path.name.removesuffix("_sample_barcode.stats.tsv")
            with path.open() as handle:
                raw_counts[sample_mod] = {row[0]: int(row[1]) for row in csv.reader(handle, delimiter="\t")}["reads"]
    overview, rna_plots, dna_plots, routed = [], [], [], {}
    for t in targets:
        mod, split = t["modality"], t["id"]
        rna = mod == "rna"
        scope = " · ".join([t["sample"], t["group"], mod.upper()] + ([] if rna else [t["mark"]]))
        path = root / f"{split}.{mod}.qc.json"
        q = json.loads(path.read_text()) if path.exists() else None
        raw = raw_counts.get(f"{t['sample']}.{mod}")
        overview.append(target_card(scope, q, rna, int(raw) if raw is not None else None))
        if q is None:
            continue
        if restart:
            routed[mod, split] = q["barcodes"] if rna else json.loads((root / f"{split}.dna.routed.json").read_text())
        if rna:
            value = saturation(q)
            if value is None:
                rna_plots.append(f'<h3>{escape(scope)}</h3><p class="unavailable">RNA saturation unavailable · no uniquely assigned GeneFull pairs.</p>')
            else:
                figure = _rna_complexity_figure(scope, [dict(branch=t["group"], value=value)])
                rna_plots.append(f'<h3>{escape(scope)}</h3>' + cards([figure], wide=True))
            note = "Region fractions use aligned blocks of uniquely mapped pairs. Exon overlap takes precedence over gene-body overlap. Strand compares R1 with the unambiguous annotation."
        else:
            dna_plots.append(f'<h3>{escape(scope)}</h3>' + dna_complexity(q, root, split, scope, t["mark"]))
            lengths = q["fragment_lengths"]
            distribution = cards([fragment_figure(lengths, scope)], wide=True)
            if any(int(length) > 1000 for length in lengths):
                distribution += details("Full fragment-length range", cards(
                    [fragment_figure(lengths, scope, full_range=True)], wide=True))
            metrics.append(details(scope + " · fragment lengths", distribution))
            note = f"Epic2 · {q['caller_status']}. Final aligned fragments are complete primary, proper pairs on canonical chromosomes after blacklist filtering, including duplicate-flagged fragments. Every eligible fragment contributes to coverage and the FRiP denominator, regardless of peak overlap. CPM coverage uses 50 bp bins."
        metrics.append(details(scope + " · analysis metrics", table(METRIC_HEADERS, rna_metrics(q) if rna else dna_metrics(q))
                               + f"<p>{escape(note)}</p>"))

    if restart:
        panels = [p for p in routed_panels(targets, routed) if p.tag != "SB"]
    library = ", ".join(sorted({t["library"] for t in targets}))
    stage = ("Measurements start at the supplied routed pairs. Raw depth, upstream losses, unassigned totals "
             "and percentages of raw input are unavailable." if restart else
             "Raw counts cover the whole sample and modality. Final counts and percentages describe each group or DNA mark. Sample barcodes use raw input pairs.")
    warnings_html = "".join(f'<li>{escape(w)}</li>' for w in warnings)
    if warnings:
        warnings_html = '<aside class="bulk-warnings"><strong>Unassigned sample barcodes</strong><ul>' + warnings_html + '</ul></aside>'
    css = """
    .run-plots{min-width:0}.run-plots>.plot-card{min-width:0}
    .bulk-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(340px,100%),1fr));gap:.8rem;margin:.75rem 0 1.2rem}
    .bulk-grid .plot-card{min-width:0;padding:.65rem}.bulk-grid .svg-scroll svg{width:100%;max-width:560px;height:auto}
    .bulk-grid .plot-actions{margin-bottom:.2rem}.bulk-warnings{border-left:3px solid var(--orange);padding:.5rem 1rem;background:#fff6ed;margin:1rem 0}
    .bulk-warnings ul{margin:.25rem 0}section{scroll-margin-top:3rem}.bulk-grid+h3{margin-top:1.6rem}
    """
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{escape(library)} · Bulk TrES-seq QC</title><style>{CSS}{css}</style></head>
<body><header class="hero"><div class="brand">TrESFlow</div><h1>Bulk TrES-seq QC</h1></header>
<nav aria-label="Report sections"><a href="#overview">Overview</a><a href="#barcodes">Barcode abundance</a>
<a href="#rna-complexity">RNA saturation</a><a href="#complexity">DNA complexity</a><a href="#metrics">Detailed metrics</a></nav><main>
<section id="overview"><h2>{escape(library)}</h2><p class="subtle">{stage}</p>{warnings_html}
<div class="run-cards">{"".join(overview)}</div></section>
<section id="barcodes"><h2>Barcode abundance</h2>{'<p class="unavailable">Raw-input sample barcode composition unavailable for supplied FASTQs.</p>' if restart else ""}{barcode_section(panels)}</section>
<section id="rna-complexity"><h2>RNA sequencing saturation</h2><p>Fraction of uniquely assigned GeneFull pairs that do not contribute a new corrected molecule · 1 − molecules / assigned pairs. UMI correction and gene assignments are fixed by STARsolo.</p>{"".join(rna_plots) if rna_plots else '<p class="unavailable">No RNA analysis.</p>'}</section>
<section id="complexity"><h2>DNA library complexity</h2><p>Picard pairs examined are split into observed unique, PCR/library duplicate, and optical duplicate pairs. Predictions use Picard’s ROI histogram and describe unique pairs expected from deeper sequencing of the same library. Duplicate-flagged fragments remain included in bulk coverage and peaks.</p>{"".join(dna_plots) if dna_plots else '<p class="unavailable">No DNA analysis.</p>'}</section>
<section id="metrics"><h2>Detailed metrics</h2>{"".join(metrics)}</section></main>{PLOT_EXPORT_SCRIPT}</body></html>'''

if __name__ == "__main__":
    spec = json.loads(base64.b64decode(sys.argv[1]))
    root = Path(sys.argv[2])
    paths = [root / f"{target['id']}.counts.tsv" for target in spec["targets"] if target["modality"] == "rna"]
    paths = [path for path in paths if path.is_file()]
    if paths:
        merge_counts(paths, Path("gene_by_group.tsv"))
    Path("bulk_report.html").write_text(render(spec, root))
