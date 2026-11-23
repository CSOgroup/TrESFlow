#!/usr/bin/env python3
"""Render pilot QC from declared process outputs and combine raw integer RNA counts."""
import base64
import csv
import json
import sys
from html import escape
from itertools import zip_longest
from contextlib import ExitStack
from pathlib import Path


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
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def metric(label, value, denominator=None, scope=""):
    pct = f"{100 * value / denominator:.2f}%" if denominator and value is not None else "—"
    return label, value, pct, scope


def plot(points, title, xlabel, ylabel):
    width, height = 440, 110
    xmax = max([p[0] for p in points] + [1])
    ymax = max([p[1] for p in points] + [1])
    coords = " ".join(f"{50 + x * width / xmax:.2f},{height + 10 - y * height / ymax:.2f}" for x, y in points)
    return (f"<figure><figcaption>{escape(title)}</figcaption><svg viewBox='0 0 510 160' role='img' "
            f"aria-label='{escape(title)}'><path d='M50 10V120H490' stroke='#aaa' fill='none'/>"
            f"<polyline points='{coords}' fill='none' stroke='#256b74' stroke-width='2'/>"
            f"<text x='0' y='12'>{ymax:,.0f}</text><text x='38' y='134'>0</text>"
            f"<text x='490' y='134' text-anchor='end'>{xmax:,.0f}</text>"
            f"<text x='250' y='152' text-anchor='middle'>{escape(xlabel)}</text></svg>"
            f"<small>{escape(ylabel)}</small></figure>")


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


def render(spec, root):
    sections = ["<h1>Bulk TrES-seq QC</h1>"]
    if spec["restart"]:
        sections.append("<p>Measurements start at the supplied routed pairs. Raw depth, upstream losses, "
                        "unassigned totals and percentages of raw input are unavailable.</p>")
    else:
        for path in sorted(root.glob("*_barcode_composition.tsv")):
            with path.open() as h:
                composition = list(csv.DictReader(h, delimiter="\t"))
            sample, modality = composition[0]["sample_id"], composition[0]["modality"]
            prefix = f"{sample}.{modality}"
            stats = (root / f"{prefix}_sample_barcode.stats.tsv").read_text().splitlines()
            raw = int(next(line.split()[1] for line in stats if line.startswith("reads\t")))
            sections.append(f"<h2>{escape(sample)} · {modality.upper()} preprocessing</h2>")
            no_sb = int(next(line.split()[1] for line in stats if line.startswith("reads_without_bc\t")))
            gates = [metric("Raw input pairs", raw),
                     metric("Unassigned at sample-barcode tagging", no_sb, raw, "raw input pairs")]
            if no_sb:
                sections.append(f"<p><b>Unassigned</b> · {no_sb:,} raw pairs had no accepted sample barcode at tagging.</p>")
            with (root / f"{prefix}_barcode_gates.tsv").open() as h:
                for row in csv.DictReader(h, delimiter="\t"):
                    gates.append(metric(row["metric"].replace("_", " "), int(row["pairs"]), raw, "raw input pairs"))
            sections.append(table(["Stage", "Read pairs", "Percent", "Denominator"], gates))
            for row in composition:
                if row["barcode_type"] == "sample_barcode" and row["category"] == "NoMatch" and int(row["count"]):
                    sections.append(f"<p><b>Unassigned</b> · {int(row['count']):,} ligation-accepted split-input pairs had no sample-barcode group.</p>")
            sections.append("<details><summary>Split-input barcode composition</summary>" +
                            table(["Barcode", "Group", "Sequence or label", "Read pairs", "Percent", "Denominator"],
                                  [(r["barcode_type"], r["group"], r["category"], int(r["count"]),
                                    f"{float(r['percentage']):.2f}%", r["denominator_definition"])
                                   for r in composition]) + "</details>")
    for target in sorted(spec["targets"], key=lambda t: (t["sample"], t["group"], t["modality"], t["id"])):
        mod, split = target["modality"], target["id"]
        path = root / f"{split}.{mod}.qc.json"
        title = " · ".join(str(t) for t in [target["library"], target["sample"], target["group"],
                                           mod.upper(), target["mark"] if mod == "dna" else None] if t)
        sections.append(f"<h2>{escape(title)}</h2>")
        if not path.exists():
            sections.append("<p>No routed analysis output was observed.</p>")
            continue
        q = json.loads(path.read_text())
        rna = mod == "rna"
        sections.append(table(["Metric", "Count", "Percent", "Denominator"], rna_metrics(q) if rna else dna_metrics(q)))
        if spec["restart"]:
            barcodes = q["barcodes"] if rna else json.loads((root / f"{split}.dna.routed.json").read_text())
            rows = []
            for tag, values in sorted(barcodes.items()):
                n = sum(values.values())
                rows.extend(metric(f"{tag} · {seq}", v, n, "supplied routed pairs") for seq, v in sorted(values.items()))
            sections.append("<details><summary>Barcode observations on routed pairs</summary>" +
                            table(["Barcode", "Read pairs", "Percent", "Denominator"], rows) + "</details>")
        if rna:
            sections.append("<p>Region fractions use aligned blocks of uniquely mapped pairs. Exon overlap takes "
                            "precedence over gene-body overlap. Strand compares R1 with the unambiguous annotation.</p>")
            labels = ["corrected molecules", "genes"]
        else:
            sections.append(f"<p>Epic2 · {escape(q['caller_status'])}. CPM coverage in 50 bp bins includes duplicate-flagged target fragments.</p>")
            # Bin fragment lengths for display only. The JSON retains exact lengths.
            bins = {}
            for length, count in q["fragment_lengths"].items():
                lo = int(length) // 25 * 25
                bins[lo] = bins.get(lo, 0) + count
            sections.append(plot(sorted(bins.items()), "Fragment lengths (25 bp bins)", "Fragment length (bp)", "Y axis · fragments"))
            labels = ["usable fragments"]
        for i, label in enumerate(labels, 1):
            sections.append(plot([(0, 0)] + [(v[0], v[i]) for v in q["depth"]], f"Depth and {label}",
                                 "Retained sequencing pairs" if rna else "Sequencing pairs at BWA", "Y axis · " + label))
        sections.append("<p>Expected yield under independent thinning of observed pairs at 10, 25, 50, 75 and 100%. "
                        + ("Uses observed support per corrected molecule and per gene with assignments fixed." if rna else
                           "Fragment yield scales linearly with retained depth under the target duplicate policy.") + "</p>")
    return ("<!doctype html><html lang='en'><meta charset='utf-8'><title>Bulk TrES-seq QC</title>"
            "<style>body{font:16px system-ui;color:#24363b;max-width:64em;margin:3em auto;padding:0 1em}"
            "table{border-collapse:collapse}th,td{padding:.3em .8em;border-bottom:1px solid #ddd;text-align:left}"
            "h2{margin-top:2em}details{margin:1em 0}figure{display:inline-block;width:27em;max-width:95%;margin:1em}"
            "svg{width:100%;font:11px system-ui}figcaption{font-weight:bold}</style>" + "\n".join(sections) + "</html>")


if __name__ == "__main__":
    spec = json.loads(base64.b64decode(sys.argv[1]))
    root = Path(sys.argv[2])
    paths = list(root.glob("*.counts.tsv"))
    if paths:
        merge_counts(paths, Path("gene_by_group.tsv"))
    Path("bulk_report.html").write_text(render(spec, root))
