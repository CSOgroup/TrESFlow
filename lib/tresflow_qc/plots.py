"""Deterministic, offline SVG figures generated from the normalized QC model.

The report embeds these SVG strings directly.  The standalone assessor can
write the exact same strings as useful standalone SVG files, so there is one
visualization implementation and no report-process plot artifact.
"""

from __future__ import annotations

import hashlib
import math
import re
import unicodedata
from collections import defaultdict
from html import escape
from pathlib import Path
from typing import Iterable

from .core import ReportModel, numeric, retention_display_stage


BACKGROUND = "#fbfaf7"
CARD = "#ffffff"
INK = "#243238"
MUTED = "#697577"
GRID = "#dfe5e2"
RNA = "#248a82"
DNA = "#d77a2b"
NOMATCH = "#aeb4b4"
UNIQUE = "#5ba89d"
PCR_DUP = "#dc9657"
OPTICAL_DUP = "#9a75ae"


def safe_name(value: object) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return cleaned or "run"


def count_label(value: object) -> str:
    number = int(float(value))
    if abs(number) >= 1_000_000:
        return f"{number / 1_000_000:.2f}M"
    if abs(number) >= 1_000:
        return f"{number / 1_000:.1f}k"
    return f"{number:,}"


def pct_label(value: object) -> str:
    return "n/a" if value is None else f"{float(value):.1f}%"


def _text_width(value: object, font_size: float = 14, *, weight: int = 400) -> float:
    """Estimate rendered system-font width for deterministic server-side SVG layout.

    Browser font metrics are unavailable while rendering the offline report, so
    layout decisions use Unicode-aware glyph classes rather than character-count
    or category-specific percentage thresholds.  The estimate is deliberately
    conservative and is exposed through ``data-layout-box`` attributes for tests.
    """
    units = 0.0
    for character in str(value):
        if character.isspace():
            units += 0.34
        elif unicodedata.east_asian_width(character) in {"W", "F"}:
            units += 1.0
        elif character in "ilI|!.,:;'`":
            units += 0.30
        elif character in "MW@#%&":
            units += 0.88
        elif character.isupper():
            units += 0.66
        elif character.isdigit():
            units += 0.58
        else:
            units += 0.54
    return units * font_size * (1.04 if weight >= 600 else 1.0)


def _ellipsize_pixels(value: object, maximum_width: float, font_size: float, *, weight: int = 400) -> str:
    text = str(value)
    if maximum_width <= _text_width("…", font_size, weight=weight):
        return ""
    if _text_width(text, font_size, weight=weight) <= maximum_width:
        return text
    low, high = 0, len(text)
    while low < high:
        middle = (low + high + 1) // 2
        candidate = text[:middle].rstrip() + "…"
        if _text_width(candidate, font_size, weight=weight) <= maximum_width:
            low = middle
        else:
            high = middle - 1
    return text[:low].rstrip() + "…"


def _layout_attributes(x: float, y: float, width: float, height: float, scope: str) -> str:
    return (
        f'data-layout-box="{x:.2f},{y:.2f},{max(0.0, width):.2f},{max(0.0, height):.2f}" '
        f'data-layout-scope="{escape(scope, quote=True)}"'
    )


def _hue_distance(left: float, right: float) -> float:
    delta = abs(left - right) % 360.0
    return min(delta, 360.0 - delta)


def categorical_palette(labels: Iterable[str], namespace: str) -> dict[str, str]:
    """Return a deterministic, scalable palette with adjacent-label contrast."""
    all_labels = set(str(label) for label in labels)
    ordered = sorted(label for label in all_labels if label != "NoMatch")
    target_separation = min(34.0, max(8.0, 300.0 / max(1, len(ordered))))
    assigned_hues: list[float] = []
    result: dict[str, str] = {}
    for label in ordered:
        digest = hashlib.sha256(f"{namespace}:{label}".encode("utf-8")).digest()
        seed = int.from_bytes(digest[:2], "big") % 360
        candidates = [(seed + 137.507764 * attempt) % 360 for attempt in range(24)]
        hue = max(
            candidates,
            key=lambda candidate: min((_hue_distance(candidate, existing) for existing in assigned_hues), default=360.0),
        )
        # Prefer the identity-derived seed whenever it is already sufficiently
        # distinct; otherwise use the deterministic maximin candidate.
        if min((_hue_distance(seed, existing) for existing in assigned_hues), default=360.0) >= target_separation:
            hue = float(seed)
        assigned_hues.append(hue)
        saturation = 46 + digest[2] % 13
        lightness = 58 + digest[3] % 11
        result[label] = f"hsl({hue:.1f} {saturation}% {lightness}%)"
    if "NoMatch" in all_labels:
        result["NoMatch"] = NOMATCH
    return result


def mark_palette(mark_names: Iterable[str]) -> dict[str, str]:
    return categorical_palette(mark_names, "mark")


def _short(value: object, length: int = 30) -> str:
    text = str(value)
    return text if len(text) <= length else text[: length - 1] + "…"


def _short_stage(stage: str) -> str:
    stage = retention_display_stage(stage)
    return {
        "After paired trimming": "Paired trimming",
        "After DT artifact filter": "DT filter",
        "L1–L2–L3 accepted": "L1–L2–L3",
        "Sample barcode accepted — all required barcodes accepted": "Sample barcode",
        "Sample barcode accepted": "Sample barcode",
        "MO barcode accepted — all required barcodes accepted": "MO barcode",
        "Mark-specific routing branch": "Mark routing",
        "Canonical marked-duplicate BAM": "Canonical markeddup",
        "Canonical NoDup final": "NoDup final",
        "Canonical chromosomes": "Canonical",
        "Called-cell final BAM": "Called-cell BAM",
        "Joint barcode accepted": "Joint barcode",
    }.get(stage, stage)


def _svg(width: int, height: int, body: str, title: str) -> str:
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(title, quote=True)}">
<title>{escape(title)}</title>
<style>
text{{font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;fill:{INK}}}
.muted{{fill:{MUTED}}}.axis{{font-size:13px}}.stage{{font-size:14px;font-weight:650}}.value{{font-size:13px;font-weight:650}}.loss{{font-size:12px;fill:{MUTED}}}.panel-title{{font-size:16px;font-weight:750}}.legend{{font-size:13px}}
</style><rect data-export-background="white" x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>{body}</svg>'''


def _retention_panel(rows: list[dict[str, object]], color: str, title: str, width: int, y_offset: int) -> str:
    rows = sorted(rows, key=lambda row: int(row["stage_order"]))
    if not rows:
        return ""
    panel_height = 272
    values = [int(row["output_pairs"]) for row in rows]
    maximum = max(values + [1])
    y_tick_width = max(_text_width(count_label(maximum * fraction), 13) for fraction in (0.0, 0.5, 1.0))
    maximum_stage_width = max((_text_width(_short_stage(str(row["stage"])), 14, weight=650) for row in rows), default=0.0)
    maximum_value_width = max((_text_width(count_label(value), 13, weight=650) for value in values), default=0.0)
    side_margin = max(34.0, maximum_stage_width / 2 + 20.0, maximum_value_width / 2 + 20.0)
    left, right = max(112.0, y_tick_width + 54.0, side_margin), width - side_margin
    top, bottom = y_offset + 78.0, y_offset + 177.0
    count = len(rows)
    xs = [left + (right - left) * index / max(1, count - 1) for index in range(count)]
    ys = [bottom - (bottom - top) * value / maximum for value in values]
    parts = [
        f'<rect x="8" y="{y_offset + 4}" width="{width - 16}" height="{panel_height - 10}" rx="10" fill="{CARD}" stroke="{GRID}"/>',
        f'<text class="panel-title" x="26" y="{y_offset + 31}" {_layout_attributes(26, y_offset + 13, _text_width(title, 17, weight=750), 22, f"retention-title-{y_offset}")}>{escape(title)}</text>',
    ]
    for fraction_value in (0.0, 0.5, 1.0):
        y = bottom - (bottom - top) * fraction_value
        label = count_label(round(maximum * fraction_value))
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="{GRID}" stroke-width="1"/>')
        parts.append(f'<text class="axis muted" x="{left - 10}" y="{y + 5:.1f}" text-anchor="end">{escape(label)}</text>')
    parts.append(f'<text class="axis muted" transform="translate(20 {y_offset + 125}) rotate(-90)" text-anchor="middle">Surviving pairs</text>')
    path = [f"M {xs[0]:.1f} {ys[0]:.1f}"]
    for index in range(1, count):
        path.extend([f"H {xs[index]:.1f}", f"V {ys[index]:.1f}"])
    parts.append(f'<path d="{" ".join(path)}" fill="none" stroke="{color}" stroke-width="3"/>')
    raw = values[0]
    node_spacing = (right - left) / max(1, count - 1)
    for index, (row, x, y) in enumerate(zip(rows, xs, ys)):
        cumulative = float(row.get("cumulative_raw_pct") or (values[index] / raw * 100 if raw else 0))
        loss_detail = ""
        if index:
            lost = values[index - 1] - values[index]
            lost_pct = lost / values[index - 1] * 100 if values[index - 1] else 0.0
            loss_detail = f"; loss from preceding stage {lost:,} pairs ({lost_pct:.1f}%)"
        internal_stage = str(row["stage"])
        display_stage = retention_display_stage(internal_stage)
        counter_detail = (
            f'; underlying accounting stage "{internal_stage}"'
            if display_stage != internal_stage
            else ""
        )
        tooltip = (
            f"{display_stage}: {values[index]:,} pairs; {cumulative:.1f}% of raw{loss_detail}; "
            f"source {row.get('count_source', 'unavailable')}{counter_detail}"
        )
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{CARD}" stroke="{color}" stroke-width="3">'
            f'<title>{escape(tooltip)}</title></circle>'
        )
        value_y = max(y_offset + 52, y - 28)
        value_text = count_label(values[index])
        percent_text = f"{cumulative:.1f}%"
        value_width = max(_text_width(value_text, 13, weight=650), _text_width(percent_text, 13))
        value_scope = f"retention-values-{y_offset}"
        parts.append(f'<g {_layout_attributes(x - value_width / 2, value_y - 14, value_width, 32, value_scope)}>')
        parts.append(f'<text class="value" x="{x:.1f}" y="{value_y:.1f}" text-anchor="middle">{escape(value_text)}</text>')
        parts.append(f'<text class="axis muted" x="{x:.1f}" y="{value_y + 16:.1f}" text-anchor="middle">{percent_text}</text></g>')
        full_stage = _short_stage(str(row["stage"]))
        display_stage = _ellipsize_pixels(full_stage, max(72.0, node_spacing - 18.0), 14, weight=650)
        stage_width = _text_width(display_stage, 14, weight=650)
        parts.append(
            f'<text class="stage" x="{x:.1f}" y="{bottom + 32}" text-anchor="middle" '
            f'{_layout_attributes(x - stage_width / 2, bottom + 17, stage_width, 20, f"retention-stages-{y_offset}")}>{escape(display_stage)}'
            f'<title>{escape(full_stage)}</title></text>'
        )
        # Transition losses remain in the node tooltip.  This deliberately
        # prioritizes the exact stage count and cumulative percentage when
        # vertical or horizontal room is insufficient for a second annotation.
    return "".join(parts)


def _retention_figures(model: ReportModel) -> list[dict[str, object]]:
    by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in model.retention:
        by_run[str(row["run"])].append(row)
    colors = mark_palette(str(row["branch"]) for row in model.retention if row["modality"] == "DNA")
    figures: list[dict[str, object]] = []
    for run, run_rows in sorted(by_run.items()):
        rna = [row for row in run_rows if row["modality"] == "RNA"]
        dna = [row for row in run_rows if row["modality"] == "DNA"]
        branches = sorted({str(row["branch"]) for row in dna})
        panels: list[tuple[list[dict[str, object]], str, str]] = []
        if rna:
            rna_branch = sorted({str(row["branch"]) for row in rna})[0]
            panels.append(([row for row in rna if row["branch"] == rna_branch], RNA, "RNA"))
        if dna:
            template = [row for row in dna if row["branch"] == branches[0]]
            shared = [row for row in template if row["stage_scope"] == "shared"]
            panels.append((shared, DNA, "DNA shared processing"))
            for branch in branches:
                branch_rows = [row for row in dna if row["branch"] == branch and row["stage_scope"] == "branch"]
                if shared and branch_rows:
                    anchor = dict(shared[-1])
                    anchor["stage_order"] = 0
                    branch_rows = [anchor, *branch_rows]
                panels.append((branch_rows, colors[branch], f"DNA · {branch}"))
        maximum_stages = max((len(rows) for rows, _, _ in panels), default=1)
        maximum_stage_width = max(
            (_text_width(_short_stage(str(row["stage"])), 14, weight=650) for rows, _, _ in panels for row in rows),
            default=90.0,
        )
        maximum_value_width = max(
            (_text_width(count_label(row["output_pairs"]), 13, weight=650) for rows, _, _ in panels for row in rows),
            default=60.0,
        )
        step_spacing = max(128.0, maximum_stage_width + 24.0, maximum_value_width + 24.0)
        left_margin = max(112.0, maximum_stage_width / 2 + 20.0, maximum_value_width / 2 + 20.0)
        right_margin = max(34.0, maximum_stage_width / 2 + 20.0, maximum_value_width / 2 + 20.0)
        content_width = left_margin + max(0, maximum_stages - 1) * step_spacing + right_margin
        title_width = max((_text_width(title, 17, weight=750) for _, _, title in panels), default=0.0) + 60.0
        width = int(math.ceil(max(820.0, content_width, title_width)))
        height = 272 * len(panels) + 8
        body = "".join(
            _retention_panel(rows, color, title, width, index * 272)
            for index, (rows, color, title) in enumerate(panels)
        )
        figures.append(
            {
                "kind": "retention",
                "run": run,
                "title": f"{run} read-pair retention",
                "svg": _svg(width, height, body, f"{run} read-pair retention"),
            }
        )
    return figures


def _composition_label(row: dict[str, object]) -> str:
    if row["category"] == "NoMatch":
        return "NoMatch"
    label = str(row.get("barcode_label") or row.get("barcode_sequence") or row["category"])
    sequence = str(row.get("barcode_sequence") or "")
    return label if not sequence or sequence == label else f"{sequence} · {label}"


def _composition_identity(row: dict[str, object]) -> str:
    if row["category"] == "NoMatch":
        return "NoMatch"
    return str(row.get("barcode_sequence") or row.get("barcode_label") or row["category"])


def _flow_legend(
    identities: list[str],
    display_labels: dict[str, str],
    colors: dict[str, str],
    *,
    width: float,
    y_start: float,
    scope: str,
) -> tuple[list[str], float]:
    """Lay out legend items using estimated rendered widths and row wrapping."""
    margin, gap, row_height = 24.0, 20.0, 25.0
    available = width - margin * 2
    x, y = margin, y_start
    parts: list[str] = []
    for identity in identities:
        full_label = display_labels[identity]
        maximum_text = min(360.0, max(130.0, available * 0.42))
        label = _ellipsize_pixels(full_label, maximum_text, 13)
        text_width = _text_width(label, 13)
        item_width = 14.0 + 7.0 + text_width
        if x > margin and x + item_width > width - margin:
            x = margin
            y += row_height
        box_width = min(item_width, available)
        parts.append(
            f'<g data-category="{escape(identity, quote=True)}" '
            f'{_layout_attributes(x, y - 3, box_width, 18, f"{scope}-row-{int((y - y_start) / row_height)}")}> '
            f'<rect x="{x}" y="{y}" width="14" height="14" rx="2" fill="{colors[identity]}"/>'
            f'<text class="legend" x="{x + 21}" y="{y + 12}">{escape(label)}<title>{escape(full_label)}</title></text></g>'
        )
        x += item_width + gap
    return parts, y + row_height


def _static_value_panel(
    grouped: list[tuple[str, list[dict[str, object]]]],
    identities: list[str],
    display_labels: dict[str, str],
    colors: dict[str, str],
    *,
    width: float,
    y_start: float,
    scope: str,
) -> tuple[list[str], float]:
    """Render a measured category-value panel for every configured category.

    This is intentionally more than a legend: each visible row carries its
    exact count and percentage for every available modality.  It provides a
    static identification path for narrow and adjacent stacked-bar segments.
    """
    category_values = []
    for identity in identities:
        values = []
        for label, rows in grouped:
            row = next((item for item in rows if _composition_identity(item) == identity), None)
            if row is None:
                continue
            count = int(row["count"])
            denominator = int(row["denominator_count"])
            values.append(f"{label} {count:,} ({count / denominator * 100 if denominator else 0:.2f}%)")
        category_values.append((identity, values))
    if not category_values:
        return [], y_start

    margin, gap, row_height = 24.0, 18.0, 28.0
    value_texts = {}
    for identity, values in category_values:
        display = display_labels[identity]
        category = display if display == identity or display.startswith(f"{identity} ·") else f"{identity} · {display}"
        value_texts[identity] = f"{category} — {' · '.join(values)}"
    natural_width = max((_text_width(text, 13) + 42.0 for text in value_texts.values()), default=280.0)
    item_width = min(max(300.0, natural_width), max(300.0, width - 2 * margin))
    columns = max(1, int((width - 2 * margin + gap) // (item_width + gap)))
    item_width = (width - 2 * margin - gap * (columns - 1)) / columns
    parts = [f'<text class="value" x="{margin}" y="{y_start + 13}">Category values</text>']
    first_row_y = y_start + 28.0
    for index, (identity, _values) in enumerate(category_values):
        column = index % columns
        row_index = index // columns
        x = margin + column * (item_width + gap)
        y = first_row_y + row_index * row_height
        available = item_width - 28.0
        display = _ellipsize_pixels(value_texts[identity], available, 13)
        parts.append(
            f'<g data-static-category="{escape(identity, quote=True)}" '
            f'{_layout_attributes(x, y - 4, item_width, 21, f"{scope}-row-{row_index}")}>'
            f'<rect x="{x}" y="{y}" width="14" height="14" rx="3" fill="{colors[identity]}"/>'
            f'<text class="legend" x="{x + 21}" y="{y + 12}">{escape(display)}'
            f'<title>{escape(value_texts[identity])}</title></text></g>'
        )
    rows = math.ceil(len(category_values) / columns)
    return parts, first_row_y + rows * row_height


def _stacked_composition_svg(
    grouped: list[tuple[str, list[dict[str, object]]]],
    *,
    title: str,
    denominator_label: str,
    exclude_nomatch: bool,
    color_namespace: str,
    show_panel_title: bool = True,
    color_map: dict[str, str] | None = None,
    static_values: bool = False,
    internal_segment_labels: bool = True,
) -> str:
    visible_rows = [
        row
        for _, rows in grouped
        for row in rows
        if not (exclude_nomatch and row["category"] == "NoMatch")
    ]
    identities = sorted(
        {_composition_identity(row) for row in visible_rows},
        key=lambda identity: (identity == "NoMatch", identity),
    )
    display_labels: dict[str, str] = {}
    for row in visible_rows:
        display_labels.setdefault(_composition_identity(row), _composition_label(row))
    color_map = color_map or categorical_palette(identities, color_namespace)
    run_label_width = max((_text_width(label, 15, weight=650) for label, _ in grouped), default=100.0)
    run_label_area = min(390.0, max(150.0, run_label_width + 22.0))
    left = 24.0 + run_label_area
    # Larger category sets receive more horizontal space, while labels for tiny
    # segments move to the measured flow legend instead of overlapping the bar.
    bar_width = max(620.0, min(1480.0, 520.0 + 34.0 * len(identities)))
    width = int(math.ceil(left + bar_width + 34.0))
    bar_height, row_gap = 46.0, 84.0
    rows_top = 52.0 if show_panel_title else 18.0
    legend_y = rows_top + row_gap * len(grouped) + 15.0
    if static_values:
        legend_parts, legend_bottom = _static_value_panel(
            grouped,
            identities,
            display_labels,
            color_map,
            width=width,
            y_start=legend_y,
            scope=f"values-{safe_name(color_namespace)}",
        )
    else:
        legend_parts, legend_bottom = _flow_legend(
            identities,
            display_labels,
            color_map,
            width=width,
            y_start=legend_y,
            scope=f"legend-{safe_name(color_namespace)}",
        )
    height = int(math.ceil(legend_bottom + 10.0))
    parts = [f'<text class="panel-title" x="20" y="28">{escape(title)}</text>'] if show_panel_title else []
    for row_index, (label, rows) in enumerate(grouped):
        visible = [row for row in rows if not (exclude_nomatch and row["category"] == "NoMatch")]
        denominator = sum(int(row["count"]) for row in visible) if exclude_nomatch else int(rows[0]["denominator_count"])
        y = rows_top + row_index * row_gap
        shown_label = _ellipsize_pixels(label, run_label_area - 18.0, 15, weight=650)
        shown_width = _text_width(shown_label, 15, weight=650)
        parts.append(
            f'<text class="stage" x="{left - 14}" y="{y + 28}" text-anchor="end" '
            f'{_layout_attributes(left - 14 - shown_width, y + 11, shown_width, 20, f"composition-run-labels-{safe_name(color_namespace)}")}>{escape(shown_label)}'
            f'<title>{escape(label)}</title></text>'
        )
        x = left
        identity_order = {identity: index for index, identity in enumerate(identities)}
        for row in sorted(visible, key=lambda item: identity_order[_composition_identity(item)]):
            count = int(row["count"])
            value = count / denominator * 100 if denominator else 0.0
            segment = bar_width * value / 100.0
            category = _composition_label(row)
            identity = _composition_identity(row)
            color = color_map[identity]
            tooltip = f"{label}; {category}: {count:,} pairs ({value:.2f}% of {denominator_label}; total {denominator:,})"
            parts.append(
                f'<rect data-category="{escape(identity, quote=True)}" data-segment-bar="true" x="{x:.2f}" y="{y}" width="{segment:.2f}" height="{bar_height}" fill="{color}">'
                f'<title>{escape(tooltip)}</title></rect>'
            )
            value_text = f"{value:.1f}% · {count_label(count)}"
            available_width = max(0.0, segment - 12.0)
            category_text = _ellipsize_pixels(category, available_width, 13, weight=700)
            required_width = max(_text_width(category_text, 13, weight=700), _text_width(value_text, 12)) + 12.0
            # This is intentionally pixel-driven: no barcode identity or fixed
            # percentage threshold decides whether a segment receives text.
            if internal_segment_labels and count > 0 and category_text and segment >= required_width:
                text_color = INK
                label_width = max(_text_width(category_text, 13, weight=700), _text_width(value_text, 12))
                parts.append(
                    f'<g data-category="{escape(identity, quote=True)}" data-in-bar-label="true" '
                    f'{_layout_attributes(x + segment / 2 - label_width / 2, y + 7, label_width, 34, f"composition-labels-{safe_name(color_namespace)}-{row_index}")}>'
                    f'<text x="{x + segment / 2:.1f}" y="{y + 19}" text-anchor="middle" style="font-size:13px;font-weight:700;fill:{text_color}">{escape(category_text)}</text>'
                    f'<text x="{x + segment / 2:.1f}" y="{y + 36}" text-anchor="middle" style="font-size:12px;fill:{text_color}">{escape(value_text)}</text></g>'
                )
            x += segment
        parts.append(f'<rect x="{left}" y="{y}" width="{bar_width}" height="{bar_height}" fill="none" stroke="{GRID}"/>')
        denominator_text = f"Total: {denominator:,} pairs"
        parts.append(f'<text class="axis muted" x="{left}" y="{y + bar_height + 18}">{escape(denominator_text)}</text>')
    parts.extend(legend_parts)
    return _svg(width, height, "".join(parts), title)


def _sample_composition_figures(model: ReportModel) -> list[dict[str, object]]:
    exact_rows = [
        row
        for row in model.barcode_composition
        if row["barcode_type"] == "sample_barcode"
        and str(row.get("contract_level", "")).startswith("exact_sequence")
    ]
    sample_identities = {_composition_identity(row) for row in exact_rows}
    sample_colors = categorical_palette(sample_identities, "sample-barcode")
    by_parent_modality: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in exact_rows:
        by_parent_modality[(str(row["parent_sample"]), str(row["modality"]))].append(row)
    run_parents: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for branch in model.branches:
        run_parents[branch.run][branch.modality].add(branch.parent)
    figures: list[dict[str, object]] = []
    for run, modality_parents in sorted(run_parents.items()):
        grouped: list[tuple[str, list[dict[str, object]]]] = []
        for modality in ("RNA", "DNA"):
            parents = sorted(modality_parents.get(modality, set()))
            for parent in parents:
                rows = by_parent_modality.get((parent, modality), [])
                if rows:
                    label = modality if len(parents) == 1 else f"{modality} · {parent}"
                    grouped.append((label, rows))
        if not grouped:
            continue
        title = f"{run} sample barcode composition"
        figures.append(
            {
                "kind": "sample_composition",
                "run": run,
                "title": title,
                "svg": _stacked_composition_svg(
                    grouped,
                    title=title,
                    denominator_label="pairs evaluated by sample-barcode tagging",
                    exclude_nomatch=False,
                    color_namespace=f"sample-barcode-{safe_name(run)}",
                    color_map=sample_colors,
                    static_values=True,
                    internal_segment_labels=False,
                    show_panel_title=False,
                ),
            }
        )
    return figures


def _mark_composition_figures(model: ReportModel) -> list[dict[str, object]]:
    by_run_group: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in model.barcode_composition:
        if row["barcode_type"] == "dna_mark":
            by_run_group[(str(row["run"]), str(row["group"]))].append(row)
    matched_marks = {
        _composition_identity(row)
        for rows in by_run_group.values()
        for row in rows
        if row["category"] != "NoMatch"
    }
    colors = categorical_palette(matched_marks, "mark")
    figures: list[dict[str, object]] = []
    runs = sorted({run for run, _group in by_run_group})
    for run in runs:
        run_groups = [(group, rows) for (row_run, group), rows in sorted(by_run_group.items()) if row_run == run]
        grouped = [("Matched marks" if len(run_groups) == 1 else group, rows) for group, rows in run_groups]
        title = f"{run} DNA mark composition"
        figures.append(
            {
                "kind": "mark_composition",
                "run": run,
                "title": title,
                "svg": _stacked_composition_svg(
                    grouped,
                    title=title,
                    denominator_label="mark-assigned DNA pairs",
                    exclude_nomatch=True,
                    color_namespace=f"mark-{safe_name(run)}",
                    color_map=colors,
                    show_panel_title=False,
                ),
            }
        )
    return figures


def _representative_ticks(values: Iterable[float], maximum: int = 7) -> list[float]:
    ordered = sorted(set(float(value) for value in values))
    if len(ordered) <= maximum:
        return ordered
    positions = {round(index * (len(ordered) - 1) / (maximum - 1)) for index in range(maximum)}
    return [ordered[index] for index in sorted(positions)]


def _complexity_figure(
    run: str,
    rows: list[dict[str, object]],
    mark_colors: dict[str, str],
) -> dict[str, object]:
    first: dict[str, dict[str, object]] = {}
    roi: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        branch = str(row["branch"])
        first.setdefault(branch, row)
        if numeric(row.get("roi_depth_multiplier")) is not None:
            roi[branch].append(row)
    branches = sorted(first)
    has_roi = bool(roi)
    roi_branches = [branch for branch in branches if roi.get(branch)]
    branch_label_width = max((_text_width(branch, 16, weight=650) for branch in branches), default=80.0)
    branch_area = min(360.0, max(120.0, branch_label_width + 18.0))
    bar_left, bar_width = 24.0 + branch_area, max(360.0, 48.0 * max(8, len(branches)))
    metric_width = max(
        (_text_width(f"PCR/library {int(row['pcr_or_library_duplicate_pairs']):,} (100.00%)", 14) for row in first.values()),
        default=230.0,
    )
    left_panel_width = bar_left + bar_width + metric_width + 55.0
    roi_depths = sorted(
        {
            float(point["roi_depth_multiplier"])
            for points in roi.values()
            for point in points
            if numeric(point.get("roi_depth_multiplier")) is not None
        }
    )
    relevant_roi_depths = [depth for depth in roi_depths if depth <= 10] or roi_depths
    roi_width = max(480.0, 72.0 * max(6, len(_representative_ticks(relevant_roi_depths)))) if has_roi else 0.0
    width = int(math.ceil(left_panel_width + (roi_width + 72.0 if has_roi else 30.0)))
    height = int(
        max(310.0, 125.0 + 88.0 * len(branches))
        + (22.0 * len(roi_branches) if has_roi else 0.0)
    )
    parts = [f'<text class="panel-title" x="20" y="28">Observed duplicate composition</text>']
    for index, branch in enumerate(branches):
        row = first[branch]
        examined = int(row["read_pairs_examined"])
        components = [
            ("Unique", int(row["observed_unique_pairs"]), UNIQUE),
            ("PCR/library", int(row["pcr_or_library_duplicate_pairs"]), PCR_DUP),
            ("Optical", int(row["optical_duplicate_pairs"]), OPTICAL_DUP),
        ]
        y = 58 + index * 86
        shown_branch = _ellipsize_pixels(branch, branch_area - 16.0, 16, weight=650)
        parts.append(f'<text class="stage" x="{bar_left - 12}" y="{y + 24}" text-anchor="end">{escape(shown_branch)}<title>{escape(branch)}</title></text>')
        x = bar_left
        for label, count, color in components:
            value = count / examined * 100 if examined else 0.0
            segment = bar_width * value / 100.0
            parts.append(f'<rect x="{x:.2f}" y="{y}" width="{segment:.2f}" height="36" fill="{color}"><title>{escape(label)}: {count:,} ({value:.3f}%)</title></rect>')
            value_text = f"{value:.1f}%"
            if count > 0 and segment >= _text_width(value_text, 13, weight=700) + 12.0:
                parts.append(f'<text x="{x + segment / 2:.1f}" y="{y + 23}" text-anchor="middle" style="font-size:13px;font-weight:700;fill:#fff">{value_text}</text>')
            x += segment
        pcr = components[1][1]
        optical = components[2][1]
        estimate = row.get("estimated_library_size")
        pcr_pct = pcr / examined * 100 if examined else 0.0
        optical_pct = optical / examined * 100 if examined else 0.0
        parts.append(f'<rect x="{bar_left}" y="{y}" width="{bar_width}" height="36" fill="none" stroke="{GRID}"/>')
        parts.append(f'<text class="axis" x="{bar_left + bar_width + 14}" y="{y + 13}">PCR/library {pcr:,} ({pcr_pct:.2f}%)</text>')
        parts.append(f'<text class="axis" x="{bar_left + bar_width + 14}" y="{y + 31}">Optical {optical:,} ({optical_pct:.2f}%)</text>')
        parts.append(f'<text class="axis muted" x="{bar_left}" y="{y + 55}">Examined {examined:,} · Picard estimated library size {"unavailable" if estimate is None else f"{int(estimate):,}"}</text>')
    legend_y = 68 + 86 * len(branches)
    for index, (label, color) in enumerate((("Observed unique", UNIQUE), ("PCR/library duplicate", PCR_DUP), ("Optical duplicate", OPTICAL_DUP))):
        x = 28 + index * 205
        parts.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{color}"/><text class="legend" x="{x + 21}" y="{legend_y + 12}">{label}</text>')
    if has_roi:
        plot_left, plot_right = left_panel_width + 48.0, width - 34.0
        roi_legend_height = 22.0 * len(roi_branches)
        plot_top, plot_bottom = 66.0 + roi_legend_height, height - 58
        all_points = [point for values in roi.values() for point in values if float(point["roi_depth_multiplier"]) <= 10]
        if not all_points:
            all_points = [point for values in roi.values() for point in values]
        max_depth = max(float(point["roi_depth_multiplier"]) for point in all_points)
        min_depth = min(float(point["roi_depth_multiplier"]) for point in all_points)
        max_y = max(float(point["roi_estimated_unique_pairs"]) for point in all_points) or 1.0
        log_min = math.log10(max(min_depth, 1e-9))
        log_max = math.log10(max_depth) if max_depth > min_depth else log_min + 1
        parts.append(
            f'<text class="panel-title" x="{plot_left}" y="28">'
            'Predicted unique pairs with deeper sequencing</text>'
        )
        for index, branch in enumerate(roi_branches):
            legend_y_roi = 50.0 + index * 22.0
            color = mark_colors[branch]
            parts.append(
                f'<line data-roi-legend-branch="{escape(branch, quote=True)}" '
                f'x1="{plot_left}" y1="{legend_y_roi - 4:.1f}" '
                f'x2="{plot_left + 22}" y2="{legend_y_roi - 4:.1f}" '
                f'stroke="{color}" stroke-width="3"/>'
                f'<text class="legend" '
                f'data-roi-legend-label="{escape(branch, quote=True)}" '
                f'x="{plot_left + 30}" y="{legend_y_roi:.1f}">{escape(branch)}</text>'
            )
        for fraction_value in (0.0, 0.5, 1.0):
            y = plot_bottom - (plot_bottom - plot_top) * fraction_value
            parts.append(f'<line x1="{plot_left}" y1="{y:.1f}" x2="{plot_right}" y2="{y:.1f}" stroke="{GRID}"/>')
            parts.append(f'<text class="axis muted" x="{plot_left - 8}" y="{y + 5:.1f}" text-anchor="end">{escape(count_label(max_y * fraction_value))}</text>')
        for depth in _representative_ticks(float(point["roi_depth_multiplier"]) for point in all_points):
            x = plot_left + (math.log10(depth) - log_min) / (log_max - log_min) * (plot_right - plot_left)
            parts.append(f'<line x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" y2="{plot_bottom}" stroke="{GRID}"/>')
            parts.append(f'<text class="axis muted" x="{x:.1f}" y="{plot_bottom + 23}" text-anchor="middle">{depth:g}×</text>')
        for branch in branches:
            points = sorted(
                [point for point in roi.get(branch, []) if float(point["roi_depth_multiplier"]) <= 10] or roi.get(branch, []),
                key=lambda point: float(point["roi_depth_multiplier"]),
            )
            coordinates = []
            for point in points:
                depth = float(point["roi_depth_multiplier"])
                value = float(point["roi_estimated_unique_pairs"])
                x = plot_left + (math.log10(depth) - log_min) / (log_max - log_min) * (plot_right - plot_left)
                y = plot_bottom - value / max_y * (plot_bottom - plot_top)
                coordinates.append((x, y, depth, value))
            if coordinates:
                points_text = " ".join(f"{x:.1f},{y:.1f}" for x, y, _, _ in coordinates)
                parts.append(f'<polyline points="{points_text}" fill="none" stroke="{mark_colors[branch]}" stroke-width="3"/>')
                for x, y, depth, value in coordinates:
                    parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{mark_colors[branch]}"><title>{escape(branch)} {depth:g}×: {value:,.0f} estimated unique pairs</title></circle>')
        parts.append(f'<text class="axis muted" x="{(plot_left + plot_right) / 2:.1f}" y="{height - 12}" text-anchor="middle">Sequencing-depth multiplier (log scale)</text>')
    title = f"{run} DNA library complexity"
    return {"kind": "complexity", "run": run, "title": title, "svg": _svg(width, height, "".join(parts), title)}


def _complexity_figures(model: ReportModel) -> list[dict[str, object]]:
    by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in model.library_complexity:
        by_run[str(row["run"])].append(row)
    colors = mark_palette(row["branch"] for row in model.library_complexity)
    return [_complexity_figure(run, rows, colors) for run, rows in sorted(by_run.items())]


def _rna_complexity_figure(run: str, rows: list[dict[str, object]]) -> dict[str, object]:
    unique_rows: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        key = (str(row.get("branch") or "RNA"), str(row.get("source") or ""))
        unique_rows.setdefault(key, row)
    ordered = [unique_rows[key] for key in sorted(unique_rows)]
    labels = [str(row.get("branch") or "RNA") for row in ordered]
    label_width = max((_text_width(label, 14, weight=650) for label in labels), default=50.0)
    left = max(110.0, label_width + 42.0)
    bar_width = 650.0
    width = int(math.ceil(left + bar_width + 34.0))
    row_height = 76.0
    height = int(math.ceil(58.0 + len(ordered) * row_height + 42.0))
    parts = ['<text class="panel-title" x="20" y="28">RNA sequencing saturation</text>']
    for index, row in enumerate(ordered):
        value = numeric(row.get("value"))
        if value is None or value < 0 or value > 100:
            raise ValueError(f"RNA sequencing saturation outside 0–100 for {run}")
        y = 52.0 + index * row_height
        observed = float(value)
        new_fraction = 100.0 - observed
        observed_width = bar_width * observed / 100.0
        new_width = bar_width - observed_width
        label = labels[index]
        display = _ellipsize_pixels(label, left - 34.0, 14, weight=650)
        parts.append(
            f'<text class="stage" x="{left - 12}" y="{y + 25}" text-anchor="end">{escape(display)}'
            f'<title>{escape(label)}</title></text>'
        )
        parts.append(
            f'<rect x="{left}" y="{y}" width="{observed_width:.2f}" height="38" fill="{RNA}">'
            f'<title>Already-observed UMI fraction: {observed:.2f}%</title></rect>'
            f'<rect x="{left + observed_width:.2f}" y="{y}" width="{new_width:.2f}" height="38" fill="#b9ddd7">'
            f'<title>New-UMI contribution: {new_fraction:.2f}%</title></rect>'
            f'<rect x="{left}" y="{y}" width="{bar_width}" height="38" fill="none" stroke="{GRID}"/>'
        )
        if observed_width >= _text_width(f"{observed:.1f}%", 13, weight=700) + 14.0:
            parts.append(f'<text x="{left + observed_width / 2:.1f}" y="{y + 25}" text-anchor="middle" style="font-size:13px;font-weight:700;fill:#fff">{observed:.1f}%</text>')
        if new_width >= _text_width(f"{new_fraction:.1f}%", 13, weight=700) + 14.0:
            parts.append(f'<text x="{left + observed_width + new_width / 2:.1f}" y="{y + 25}" text-anchor="middle" style="font-size:13px;font-weight:700;fill:{INK}">{new_fraction:.1f}%</text>')
        parts.append(f'<text class="axis muted" x="{left}" y="{y + 57}">STARsolo reported saturation {observed:.2f}%</text>')
    legend_y = height - 27.0
    parts.append(
        f'<rect x="{left}" y="{legend_y - 11}" width="14" height="14" rx="3" fill="{RNA}"/>'
        f'<text class="legend" x="{left + 21}" y="{legend_y}">Already-observed UMI fraction</text>'
        f'<rect x="{left + 260}" y="{legend_y - 11}" width="14" height="14" rx="3" fill="#b9ddd7"/>'
        f'<text class="legend" x="{left + 281}" y="{legend_y}">New-UMI contribution</text>'
    )
    title = f"{run} RNA library complexity"
    return {"kind": "rna_complexity", "run": run, "title": title, "svg": _svg(width, height, "".join(parts), title)}


def _rna_complexity_figures(model: ReportModel) -> list[dict[str, object]]:
    by_run: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in model.qc_metrics:
        if row.get("metric") == "sequencing_saturation_pct" and numeric(row.get("value")) is not None:
            by_run[str(row["run"])].append(row)
    return [_rna_complexity_figure(run, rows) for run, rows in sorted(by_run.items())]


def generate_inline_plots(model: ReportModel) -> list[dict[str, object]]:
    """Return every available report figure as an in-memory inline SVG."""
    return [
        *_retention_figures(model),
        *_mark_composition_figures(model),
        *_sample_composition_figures(model),
        *_complexity_figures(model),
        *_rna_complexity_figures(model),
    ]


def _figure_filename(item: dict[str, object]) -> str:
    kind = str(item["kind"])
    suffix = {
        "retention": "read_retention",
        "sample_composition": "sample_barcode_composition",
        "mark_composition": "dna_mark_composition",
        "complexity": "dna_library_complexity",
        "rna_complexity": "rna_library_complexity",
    }[kind]
    return f"{safe_name(item.get('run', 'run'))}_{suffix}.svg"


def write_inline_plots(plots: list[dict[str, object]], output_dir: Path) -> list[dict[str, object]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[dict[str, object]] = []
    for item in plots:
        path = output_dir / _figure_filename(item)
        path.write_text(str(item["svg"]) + "\n", encoding="utf-8")
        outputs.append({**item, "svg_path": path})
    return outputs


def _write_kind(model: ReportModel, output_dir: Path, kind: str) -> list[dict[str, object]]:
    return write_inline_plots([item for item in generate_inline_plots(model) if item["kind"] == kind], output_dir)


def generate_retention_plots(model: ReportModel, output_dir: Path) -> list[dict[str, object]]:
    return _write_kind(model, output_dir, "retention")


def generate_sample_composition_plots(model: ReportModel, output_dir: Path) -> list[dict[str, object]]:
    return _write_kind(model, output_dir, "sample_composition")


def generate_mark_composition_plots(model: ReportModel, output_dir: Path) -> list[dict[str, object]]:
    return _write_kind(model, output_dir, "mark_composition")


def generate_complexity_plots(model: ReportModel, output_dir: Path) -> list[dict[str, object]]:
    return _write_kind(model, output_dir, "complexity")


def generate_all_plots(model: ReportModel, output_dir: Path) -> list[dict[str, object]]:
    return write_inline_plots(generate_inline_plots(model), output_dir)
