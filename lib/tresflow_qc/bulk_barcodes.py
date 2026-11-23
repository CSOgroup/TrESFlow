"""Read-pair abundance at the measurement stages available to a bulk run."""
from collections import Counter, defaultdict
from dataclasses import dataclass
from html import escape

from .plots import DNA, GRID, NOMATCH, RNA, _ellipsize_pixels, count_label, svg_figure


NAMES = {"SB": "Sample barcodes", "L1": "Ligation round L1", "L2": "Ligation round L2",
         "L3": "Ligation round L3", "MO": "DNA marks"}


@dataclass
class BarcodePanel:
    scope: str
    modality: str
    tag: str
    stage: str
    denominator: int
    definition: str
    # (sequence or mark, group or mark label) -> read pairs
    counts: Counter
    sample: str = ""

    def __post_init__(self):
        if sum(self.counts.values()) != self.denominator:
            raise ValueError(f"{self.scope} {self.tag}: barcode counts disagree with {self.definition}")

    @property
    def title(self):
        return NAMES[self.tag]


def raw_panels(sample, modality, composition, tag_counts, raw_pairs, unmatched):
    scope = f"{sample} · {modality.upper()}"
    groups = {r["category"]: r["group"] for r in composition
              if r["barcode_type"] == "sample_barcode" and r["category"] != "NoMatch"}
    tag_counts = dict(tag_counts)
    if tag_counts.setdefault("NoMatch", unmatched) != unmatched:
        raise ValueError(f"{scope}: sample-barcode NoMatch disagrees with tagging stats")
    tagged = Counter({(seq, groups.get(seq, "")): n for seq, n in tag_counts.items()})
    panels = [BarcodePanel(scope, modality, "SB", "Sample-barcode tagging", raw_pairs, "raw input pairs", tagged, sample=sample)]
    subsets = defaultdict(list)
    for row in composition:
        kind = row["barcode_type"]
        group = row["group"] if kind == "dna_mark" else ""
        subsets[kind, group].append(row)
    for (kind, group), rows in subsets.items():
        tag = {"sample_barcode": "SB", "dna_mark": "MO"}.get(kind, kind)
        counts = Counter()
        for row in rows:
            label = row["group"] if tag == "SB" and row["category"] != "NoMatch" else ""
            counts[row["category"], label] += int(row["count"])
        panels.append(BarcodePanel(scope if not group else f"{scope} · {group}", modality, tag,
                                   "Split input", int(rows[0]["denominator_count"]),
                                   rows[0]["denominator_definition"], counts, sample=sample))
    return panels


def routed_panels(targets, observations):
    """Pool disjoint routed DNA marks within each group, preserving mark names on MO."""
    grouped = defaultdict(Counter)
    totals = Counter()
    for t in targets:
        mod, split = t["modality"], t["id"]
        barcodes = observations.get((mod, split))
        if barcodes is None:
            continue
        key = (t["sample"], mod, t["group"])
        # Each supplied pair has exactly one validated corrected L1 sequence.
        totals[key] += sum(barcodes["L1"].values())
        for tag in NAMES:
            if tag == "MO" and mod != "dna":
                continue
            counts = grouped[key, tag]
            for seq, n in barcodes[tag].items():
                counts[seq, t["mark"] if tag == "MO" else ""] += n
    return [BarcodePanel(f"{sample} · {mod.upper()} · {group}", mod, tag, "Supplied routed pairs",
                         totals[sample, mod, group], "supplied routed read pairs in this group", counts, sample=sample)
            for ((sample, mod, group), tag), counts in grouped.items()]


def label(identity):
    seq, group = identity
    name = {"NoMatch": "Unassigned", "unavailable": "Tag unavailable"}.get(seq, seq)
    return f"{name} · {group}" if group else name


def abundance_svg(panel):
    special = {"NoMatch", "unavailable"}
    ranked = sorted(panel.counts.items(), key=lambda item: (-item[1], item[0]))
    assigned = [item for item in ranked if item[0][0] not in special]
    shown = assigned[:8] + [item for item in ranked if item[0][0] in special]
    width, left, right = 460, 142, 335
    body = []
    def text(x, y, value, style="axis", anchor="start"):
        body.append(f'<text class="{style}" x="{x}" y="{y}" text-anchor="{anchor}">{escape(value)}</text>')

    title = f"{panel.title} · {panel.scope} · {panel.stage}"
    text(10, 22, panel.title, "panel-title")
    text(10, 44, _ellipsize_pixels(panel.scope, width - 20, 13))
    if panel.stage == "Supplied routed pairs":
        population = "read pairs in group"
    elif panel.stage == "Sample-barcode tagging":
        population = "raw input pairs"
    elif panel.tag == "SB":
        population = "ligation-accepted pairs"
    elif panel.tag == "MO":
        population = "group-assigned pairs"
    else:
        population = "split-input pairs"
    text(10, 64, f"{panel.stage} · {panel.denominator:,} {population}", "muted axis")
    maximum = max((n for _, n in shown), default=0)
    # Keep the mark name and its observed oligo visible on short mobile labels.
    show_mark_sequence = panel.tag == "MO" and any(identity[1] for identity, _ in shown)
    row_height = 38 if show_mark_sequence else 24
    for i, (identity, n) in enumerate(shown):
        y = 85 + i * row_height
        name = label(identity)
        if show_mark_sequence:
            text(left - 7, y + 2, _ellipsize_pixels(identity[1], left - 15, 13), anchor="end")
            text(left - 7, y + 17, _ellipsize_pixels(identity[0], left - 15, 13), "muted axis", "end")
        else:
            text(left - 7, y + 4, _ellipsize_pixels(name, left - 15, 13), anchor="end")
        body.append(f'<rect x="{left}" y="{y-7}" width="{right-left}" height="12" rx="2" fill="{GRID}"/>')
        color = NOMATCH if identity[0] in special else (RNA if panel.modality == "rna" else DNA)
        length = (right - left) * n / maximum if maximum else 0
        pct = f"{100*n/panel.denominator:.1f}%" if panel.denominator else "unavailable"
        body.append(f'<rect x="{left}" y="{y-7}" width="{length:.3f}" height="12" rx="2" fill="{color}">'
                    f'<title>{escape(name)} · {n:,} read pairs · {pct}</title></rect>')
        display_count = count_label(n) if n >= 100_000 else f"{n:,}"
        text(right + 7, y + 4, f"{display_count} · {pct}", "value")
    bottom = 91 + len(shown)*row_height
    if maximum:
        text(left, bottom, "0", "muted axis")
        text((left + right) / 2, bottom, "Read pairs", "muted axis", "middle")
        text(right, bottom, count_label(maximum), "muted axis", "end")
        bottom += 22
    if not shown:
        text(10, bottom, "No read pairs observed.", "muted axis")
        bottom += 22
    if len(assigned) > 8:
        text(10, bottom, f"Top 8 of {len(assigned)} assigned categories. Full counts below.", "muted axis")
        bottom += 22
    return svg_figure(width, bottom, "".join(body), title)
