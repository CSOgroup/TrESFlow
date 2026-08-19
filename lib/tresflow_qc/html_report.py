"""Concise, self-contained HTML rendering for the TrESFlow QC model."""

from __future__ import annotations

import hashlib
import html
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

from .core import ReportModel, numeric, retention_display_stage


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def fmt_count(value: object) -> str:
    parsed = numeric(value)
    return "unavailable" if parsed is None else f"{int(round(parsed)):,}"


def fmt_pct(value: object, digits: int = 1) -> str:
    parsed = numeric(value)
    return "unavailable" if parsed is None else f"{parsed:.{digits}f}%"


def _pct(numerator: Optional[int], denominator: Optional[int]) -> Optional[float]:
    if numerator is None or denominator in {None, 0}:
        return None
    return numerator / denominator * 100.0


def _count_and_raw_pct(count: Optional[int], raw: Optional[int]) -> str:
    value = _pct(count, raw)
    if count is None or value is None:
        return "unavailable"
    return f"{count:,} ({value:.1f}% of raw)"


def _raw_final(raw: Optional[int], final: Optional[int]) -> str:
    value = _pct(final, raw)
    if raw is None or final is None or value is None:
        return "unavailable"
    return f"{raw:,} / {final:,} ({value:.1f}% of raw)"


def _table(headers: list[tuple[str, str]], rows: Iterable[dict[str, object]], classes: str = "") -> str:
    rows = list(rows)
    if not rows:
        return '<p class="unavailable">Unavailable for this run.</p>'
    head = "".join(f"<th>{esc(label)}</th>" for _, label in headers)
    body = []
    for row in rows:
        cells = []
        for key, _ in headers:
            value = row.get(key, "")
            if key == "value" and row.get("unit") == "percent":
                value = fmt_pct(value, 2)
            elif key == "value" and row.get("unit") in {"read_pairs", "cells", "fragments"}:
                value = fmt_count(value)
            elif key.endswith("_pct") or key in {"percentage", "percent_duplication"}:
                value = fmt_pct(value)
            elif key.endswith("_pairs") or key.endswith("_count") or (
                key in {"output_pairs", "estimated_library_size", "count"}
                and isinstance(value, (int, float))
            ):
                value = fmt_count(value)
            cells.append(f"<td>{esc(value)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="table-scroll"><table class="{esc(classes)}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def _run_cards(model: ReportModel) -> str:
    cards = []
    runs = sorted({branch.run for branch in model.branches})
    for run in runs:
        branches = [branch for branch in model.branches if branch.run == run]
        configured_groups = sorted({branch.group for branch in branches if branch.group})
        group_title = configured_groups[0] if len(configured_groups) == 1 else run
        modalities = sorted({branch.modality for branch in branches})
        marks = sorted({branch.branch for branch in branches if branch.modality == "DNA"})
        parents = sorted({branch.parent for branch in branches})
        rows = [row for row in model.retention if row["run"] == run]
        raw_by_modality: dict[str, Optional[int]] = {}
        for modality in modalities:
            raw_values = [
                int(row["output_pairs"])
                for row in rows
                if row["modality"] == modality and row["stage"] == "Raw input"
            ]
            raw_by_modality[modality] = raw_values[0] if raw_values else None
        final_rna = next(
            (
                int(row["output_pairs"])
                for row in rows
                if row["modality"] == "RNA" and row["stage"] == "Called-cell final BAM"
            ),
            None,
        )
        final_dna_values = [
            int(row["output_pairs"])
            for row in rows
            if row["modality"] == "DNA" and row["stage"] == "Canonical NoDup final"
        ]
        final_dna = sum(final_dna_values) if final_dna_values else None
        rna_assignment = next(
            (
                int(row["output_pairs"])
                for row in rows
                if row["modality"] == "RNA"
                and row["stage"] == "Sample barcode accepted — all required barcodes accepted"
            ),
            None,
        )
        dna_assignment = next(
            (
                int(row["output_pairs"])
                for row in rows
                if row["modality"] == "DNA"
                and row["stage"] == "MO barcode accepted — all required barcodes accepted"
            ),
            None,
        )
        complexity_first: dict[str, dict[str, object]] = {}
        for row in model.library_complexity:
            if row["run"] == run:
                complexity_first.setdefault(str(row["branch"]), row)
        examined = sum(int(row["read_pairs_examined"]) for row in complexity_first.values()) or None
        pcr_library = sum(int(row["pcr_or_library_duplicate_pairs"]) for row in complexity_first.values()) if complexity_first else None
        optical = sum(int(row["optical_duplicate_pairs"]) for row in complexity_first.values()) if complexity_first else None
        if complexity_first:
            for branch, row in complexity_first.items():
                examined_branch = int(row["read_pairs_examined"])
                pcr_branch = int(row["pcr_or_library_duplicate_pairs"])
                optical_branch = int(row["optical_duplicate_pairs"])
                unique_branch = int(row["observed_unique_pairs"])
                if min(examined_branch, pcr_branch, optical_branch, unique_branch) < 0:
                    raise ValueError(f"Negative DNA duplicate component for {run}/{branch}")
                if pcr_branch + optical_branch + unique_branch != examined_branch:
                    raise ValueError(f"DNA duplicate components do not reconcile for {run}/{branch}")
        if not examined or pcr_library is None or optical is None:
            duplicate_text = "unavailable"
        else:
            duplicate_text = (
                f"PCR/library: {pcr_library:,} ({pcr_library / examined * 100:.2f}% of examined pairs)\n"
                f"Optical: {optical:,} ({optical / examined * 100:.2f}% of examined pairs)"
            )
        saturation_rows: dict[tuple[str, str], dict[str, object]] = {}
        for row in model.qc_metrics:
            if row["run"] == run and row["metric"] == "sequencing_saturation_pct" and numeric(row.get("value")) is not None:
                saturation_rows.setdefault((str(row.get("branch") or "RNA"), str(row.get("source") or "")), row)
        if not saturation_rows:
            saturation_text = "Not available"
        elif len(saturation_rows) == 1:
            saturation = float(next(iter(saturation_rows.values()))["value"])
            if saturation < 0 or saturation > 100:
                raise ValueError(f"RNA sequencing saturation outside 0–100 for {run}")
            saturation_text = f"{saturation:.1f}%"
        else:
            saturation_text = f"{len(saturation_rows)} branch values; see RNA complexity"
        metrics = [
            ("Modalities", " + ".join(modalities)),
            ("DNA marks", ", ".join(marks) if marks else "none"),
        ]
        if "RNA" in modalities:
            metrics.extend(
                [
                    ("RNA raw / final", _raw_final(raw_by_modality.get("RNA"), final_rna)),
                    (
                        "RNA fully barcode-assigned pairs",
                        _count_and_raw_pct(rna_assignment, raw_by_modality.get("RNA")),
                        "Pairs passing L1, L2, L3 and sample-barcode gates on the same pair.",
                    ),
                    (
                        "RNA sequencing saturation",
                        saturation_text,
                        "STARsolo Summary.csv: 1 minus collapsed UMI count divided by reads assigned to a unique GeneFull feature. Multiple inputs are not averaged without their exact denominators.",
                    ),
                ]
            )
        if "DNA" in modalities:
            metrics.extend(
                [
                    ("DNA raw / final", _raw_final(raw_by_modality.get("DNA"), final_dna)),
                    (
                        "DNA fully barcode-assigned pairs",
                        _count_and_raw_pct(dna_assignment, raw_by_modality.get("DNA")),
                        "Pairs passing L1, L2, L3, sample-barcode and MO/mark-barcode gates on the same pair.",
                    ),
                    (
                        "DNA duplicate components",
                        duplicate_text,
                        "Each numerator and READ_PAIRS_EXAMINED denominator is summed across DNA marks before percentages are calculated. Optical duplicates are a subset of total duplicate pairs; PCR/library duplicates are total minus optical.",
                    ),
                ]
            )
        cells = []
        for item in metrics:
            label, value = item[:2]
            title = item[2] if len(item) > 2 else ""
            title_attribute = f' title="{esc(title)}"' if title else ""
            cells.append(f'<div{title_attribute}><span>{esc(label)}</span><strong>{esc(value)}</strong></div>')
        cards.append(
            f'''<article class="run-card">
              <h3>{esc(group_title)}</h3>
              <p class="subtle">{esc(", ".join(parents))}</p>
              <div class="mini-grid">{"".join(cells)}</div>
            </article>'''
        )
    return "".join(cards) or '<p class="unavailable">No independent-run identities could be reconstructed.</p>'


def _overview(model: ReportModel) -> str:
    modalities = sorted({branch.modality for branch in model.branches})
    groups = sorted({branch.run for branch in model.branches})
    marks = sorted({branch.branch for branch in model.branches if branch.modality == "DNA"})
    dna_modes = []
    for sample in model.run_metadata.get("samples", []):
        if sample.get("modality") == "dna":
            dna_modes.append(f"{sample.get('id')}: {sample.get('dna_tagmentation') or 'unspecified'}")
    if not dna_modes:
        dna_parents = sorted(
            {
                str(branch.parent)
                for branch in model.branches
                if branch.modality == "DNA"
            }
        )
        dt_parents = {
            str(row["parent_sample"])
            for row in model.qc_metrics
            if row["metric"] == "dual_tag_filter_input_pairs"
        }
        dna_modes = [
            f"{parent}: {'dual' if parent in dt_parents else 'unspecified'}"
            for parent in dna_parents
        ]
    return f'''<section id="overview"><h2>Run overview</h2><p class="library-name">{esc(model.library_name)}</p>
      <div class="overview-strip">
        <div><span>Independent runs</span><strong>{esc(", ".join(groups) or "unavailable")}</strong></div>
        <div><span>Modalities</span><strong>{esc(", ".join(modalities) or "unavailable")}</strong></div>
        <div><span>DNA marks</span><strong>{esc(", ".join(marks) or "none")}</strong></div>
        <div><span>DNA tagmentation</span><strong>{esc("; ".join(dna_modes) if dna_modes else "not present")}</strong></div>
        <div><span>Pipeline version</span><strong>{esc(model.pipeline_version)}</strong></div>
        <div><span>Generated</span><strong>{esc(model.generated_at)}</strong></div>
      </div>
      <div class="run-cards">{_run_cards(model)}</div>
    </section>'''


def _plots_of_kind(plots: list[dict[str, object]], kind: str) -> list[dict[str, object]]:
    return [item for item in plots if item["kind"] == kind]


def _plot_markup(item: dict[str, object]) -> str:
    return (
        f'<div class="plot-actions"><button type="button" id="{esc(item["svg_download_id"])}" '
        f'class="plot-download" data-svg-download="{esc(item["dom_id"])}" data-export-format="svg" '
        f'data-filename="{esc(item["svg_filename"])}" aria-controls="{esc(item["dom_id"])}">Download SVG</button>'
        f'<button type="button" id="{esc(item["png_download_id"])}" '
        f'class="plot-download" data-svg-download="{esc(item["dom_id"])}" data-export-format="png" '
        f'data-filename="{esc(item["png_filename"])}" aria-controls="{esc(item["dom_id"])}">Download PNG</button></div>'
        f'<div class="svg-scroll">{item["svg"]}</div>'
    )


def _retention_section(plots: list[dict[str, object]]) -> str:
    selected = _plots_of_kind(plots, "retention")
    if not selected:
        return ""
    cards = "".join(
        f'<article class="plot-card run-plot"><h3>{esc(item["run"])}</h3>{_plot_markup(item)}</article>'
        for item in selected
    )
    return f'<section id="retention"><h2>Read-pair retention</h2><div class="run-plots">{cards}</div></section>'


def _single_plot_section(plots: list[dict[str, object]], kind: str, heading: str, section_id: str) -> str:
    selected = _plots_of_kind(plots, kind)
    if not selected:
        return ""
    figures = "".join(
        f'<article class="plot-card run-plot"><h3>{esc(item["run"])}</h3>{_plot_markup(item)}</article>'
        for item in selected
    )
    return f'<section id="{esc(section_id)}"><h2>{esc(heading)}</h2><div class="run-plots">{figures}</div></section>'


def _sample_section(plots: list[dict[str, object]]) -> str:
    selected = _plots_of_kind(plots, "sample_composition")
    if not selected:
        return ""
    panels = "".join(
        f'<article class="plot-card run-plot"><h3>{esc(item["run"])}</h3>{_plot_markup(item)}</article>'
        for item in selected
    )
    return f'<section id="sample-composition"><h2>Sample barcode composition</h2><div class="run-plots">{panels}</div></section>'


def _dna_complexity_section(plots: list[dict[str, object]]) -> str:
    selected = _plots_of_kind(plots, "complexity")
    if not selected:
        return ""
    cards = "".join(
        f'<article class="plot-card run-plot"><h3>{esc(item["run"])}</h3>{_plot_markup(item)}</article>'
        for item in selected
    )
    explanation = (
        '<details class="metric-explanation"><summary>Percent duplicates</summary>'
        '<p>Picard reports total duplicate pairs among <code>READ_PAIRS_EXAMINED</code>. '
        'TrESFlow supplies coordinate-sorted BWA-MEM2 records after blacklist-overlap exclusion and a SAM proper-pair flag requirement; '
        'canonical-chromosome normalization occurs after MarkDuplicates, and TrESFlow applies no MAPQ threshold before this metric. '
        'Optical duplicates are the spatially classified subset of total duplicates. PCR/library duplicates are total duplicate pairs minus optical duplicate pairs. '
        'The denominator is the pairs Picard reports as examined after Picard’s own eligibility rules. Duplication reflects both sequencing depth and library complexity.</p></details>'
        '<details class="metric-explanation"><summary>Predicted unique pairs with deeper sequencing</summary>'
        '<p>Picard MarkDuplicates estimates the return on investment from sequencing the same library more deeply. '
        'For each sequencing-depth multiplier, Picard reports an estimated unique-yield multiplier; TrESFlow converts '
        'that value to predicted unique read pairs using the currently observed unique pairs. '
        'Curves are shown separately for each DNA mark. A curve that progressively flattens indicates diminishing '
        'recovery of new unique molecules as additional sequencing increasingly revisits molecules already observed.</p></details>'
    )
    return f'<section id="complexity"><h2>DNA library complexity</h2>{explanation}<div class="run-plots">{cards}</div></section>'


def _rna_complexity_section(plots: list[dict[str, object]]) -> str:
    selected = _plots_of_kind(plots, "rna_complexity")
    if not selected:
        return ""
    cards = "".join(
        f'<article class="plot-card run-plot"><h3>{esc(item["run"])}</h3>{_plot_markup(item)}</article>'
        for item in selected
    )
    explanation = (
        '<details class="metric-explanation"><summary>RNA sequencing saturation</summary>'
        '<p>STARsolo 2.7.11b reports <code>1 − collapsed UMIs / reads assigned to a unique GeneFull feature</code>. '
        'With TrESFlow’s STARsolo settings, UMIs are collapsed with <code>1MM_CR</code> and filtered with <code>MultiGeneUMI_CR</code>. '
        'The value is the fraction of unique-feature reads that do not contribute a new reported UMI, and depends on sequencing depth and RNA library complexity; '
        'it is not an absolute RNA duplicate-pair count.</p></details>'
    )
    return f'<section id="rna-complexity"><h2>RNA library complexity and saturation</h2>{explanation}<div class="run-plots">{cards}</div></section>'


def _retention_table(model: ReportModel) -> str:
    source_rows = sorted(
        model.retention,
        key=lambda row: (
            str(row["run"]),
            row["modality"] != "RNA",
            str(row["branch"]),
            int(row["stage_order"]),
        ),
    )
    rows = [
        {
            **row,
            "display_stage": retention_display_stage(row["stage"]),
            "accounting_stage": row["stage"],
        }
        for row in source_rows
    ]
    return _table(
        [
            ("run", "Independent run"),
            ("modality", "Modality"),
            ("branch", "Branch"),
            ("display_stage", "Stage"),
            ("accounting_stage", "Accounting stage"),
            ("output_pairs", "Surviving pairs"),
            ("retained_prev_pct", "Of preceding"),
            ("cumulative_raw_pct", "Of raw"),
            ("unit", "Unit"),
            ("count_source", "Source"),
        ],
        rows,
    )


def _composition_table(model: ReportModel) -> str:
    rows = sorted(
        model.barcode_composition,
        key=lambda row: (
            str(row["parent_sample"]),
            str(row["barcode_type"]),
            str(row["modality"]),
            str(row["group"]),
            row["category"] == "NoMatch",
            str(row["barcode_label"]),
        ),
    )
    return _table(
        [
            ("run", "Independent run"),
            ("parent_sample", "Parent sample"),
            ("modality", "Modality"),
            ("barcode_type", "Barcode type"),
            ("group", "Group"),
            ("barcode_label", "Label"),
            ("barcode_sequence", "Sequence"),
            ("count", "Pairs"),
            ("percentage", "Percent"),
            ("denominator_count", "Denominator"),
            ("denominator_definition", "Denominator definition"),
            ("source", "Source"),
        ],
        rows,
    )


def _complexity_table(model: ReportModel) -> str:
    first: dict[tuple[object, object], dict[str, object]] = {}
    for row in model.library_complexity:
        first.setdefault((row["run"], row["branch"]), row)
    return _table(
        [
            ("run", "Independent run"),
            ("branch", "Mark"),
            ("read_pairs_examined", "Examined"),
            ("observed_unique_pairs", "Observed unique"),
            ("pcr_or_library_duplicate_pairs", "PCR/library duplicates"),
            ("optical_duplicate_pairs", "Optical duplicates"),
            ("percent_duplication", "Total duplication"),
            ("estimated_library_size", "Picard estimated library size"),
        ],
        [first[key] for key in sorted(first)],
    )


def _details(model: ReportModel) -> str:
    metrics = sorted(
        model.qc_metrics,
        key=lambda row: (
            str(row["run"]),
            row["modality"] != "RNA",
            str(row["branch"]),
            str(row["metric"]),
        ),
    )
    metric_table = _table(
        [
            ("run", "Independent run"),
            ("modality", "Modality"),
            ("branch", "Branch"),
            ("metric", "Metric"),
            ("value", "Value"),
            ("unit", "Unit"),
            ("denominator", "Denominator / definition"),
            ("source", "Source"),
        ],
        metrics,
    )
    filenames = ("read_retention.tsv", "qc_metrics.tsv", "barcode_composition.tsv", "library_complexity.tsv")
    links = "".join(f'<a class="file-link" href="{name}" download>{name}</a>' for name in filenames)
    return f'''<section id="details"><h2>Detailed metrics</h2><div class="file-links">{links}</div>
      <details open><summary>Read-pair retention</summary>{_retention_table(model)}</details>
      <details><summary>Barcode composition</summary>{_composition_table(model)}</details>
      <details><summary>QC metrics</summary>{metric_table}</details>
      <details><summary>DNA library complexity</summary>{_complexity_table(model)}</details>
    </section>'''


CSS = """
:root{--paper:#f5f7f5;--card:#fff;--ink:#263433;--muted:#64716f;--line:#dce5e1;--teal:#2c8982;--teal-soft:#dcefeb;--orange:#cf782f;--orange-soft:#faeadc;--focus:#286f9b}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.52 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
header{max-width:1540px;margin:0 auto;padding:2rem max(3vw,1.2rem) 1rem}.brand{font-size:.75rem;letter-spacing:.12em;text-transform:uppercase;color:var(--teal);font-weight:760}.hero h1{font-size:clamp(1.75rem,3vw,2.65rem);line-height:1.12;margin:.4rem 0 0;color:var(--ink)}
nav{position:sticky;top:0;z-index:10;padding:.7rem max(4vw,2rem);background:rgba(245,247,245,.96);border-block:1px solid var(--line);display:flex;gap:.5rem;overflow:auto;white-space:nowrap;backdrop-filter:blur(10px)}nav a{color:var(--ink);text-decoration:none;font-size:.84rem;font-weight:680;padding:.42rem .72rem;border-radius:999px}nav a:hover,nav a:focus-visible{background:var(--teal-soft);outline:none;box-shadow:0 0 0 3px rgba(44,137,130,.18)}
main{max-width:1540px;margin:auto;padding:0 max(3vw,1.2rem) 4rem}section{padding:2.7rem 0}section+section{border-top:1px solid var(--line)}h2{font-size:1.55rem;line-height:1.2;margin:0 0 1.25rem}h3{margin:.2rem 0}.library-name{color:var(--muted);font-size:1rem;margin:-.75rem 0 1.2rem}.subtle{color:var(--muted);margin:.25rem 0}
.overview-strip{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:1px;background:var(--line);border:1px solid var(--line);border-radius:18px;overflow:hidden;margin:1rem 0 1.5rem;box-shadow:0 2px 7px rgba(35,63,58,.04)}.overview-strip div{background:var(--card);padding:1rem 1.1rem}.overview-strip strong,.overview-strip span{display:block}.overview-strip span,.mini-grid span{font-size:.69rem;color:var(--muted);text-transform:uppercase;letter-spacing:.055em}.overview-strip strong{font-size:.91rem;margin-top:.15rem;overflow-wrap:anywhere}
.run-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(390px,100%),1fr));gap:1.1rem}.run-card{background:var(--card);border:1px solid var(--line);border-radius:20px;padding:1.25rem;box-shadow:0 5px 18px rgba(35,63,58,.075)}.mini-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.85rem;margin-top:1rem}.mini-grid div{border-top:1px solid var(--line);padding-top:.6rem}.mini-grid span,.mini-grid strong{display:block}.mini-grid strong{font-size:.92rem;margin-top:.2rem;overflow-wrap:anywhere;white-space:pre-line}
.run-plots{display:grid;gap:1.5rem}.plot-card{background:var(--card);border:1px solid var(--line);border-radius:22px;padding:1.15rem;box-shadow:0 7px 22px rgba(35,63,58,.085)}.run-plot h3{font-size:1.12rem;margin:0 0 .65rem;padding-bottom:.65rem;border-bottom:1px solid var(--line)}.plot-actions{display:flex;justify-content:flex-end;flex-wrap:wrap;gap:.45rem;margin:0 0 .55rem}.plot-download{border:0;border-radius:999px;background:var(--teal-soft);color:#145f5a;font:700 .78rem/1 system-ui;padding:.58rem .8rem;cursor:pointer}.plot-download:hover{background:#cbe6e1}.plot-download:focus-visible{outline:none;box-shadow:0 0 0 3px rgba(40,111,155,.28)}.svg-scroll{width:100%;overflow-x:auto;padding:.2rem 0}.svg-scroll svg{display:block;width:auto;height:auto;max-width:none;margin-inline:auto}
.metric-explanation{background:var(--card);border:1px solid var(--line);border-radius:16px;padding:.75rem 1rem;margin:-.35rem 0 1.2rem;box-shadow:0 2px 8px rgba(35,63,58,.04)}.metric-explanation summary{cursor:pointer;font-weight:720;color:#176b66}.metric-explanation p{max-width:1000px;color:var(--muted);margin:.7rem 0 .25rem}
.table-scroll{overflow:auto;max-height:540px;border:1px solid var(--line);border-radius:14px;margin:.8rem 0}table{border-collapse:collapse;width:100%;font-size:.79rem;background:var(--card)}th{position:sticky;top:0;background:#edf3f0;text-align:left;z-index:1}th,td{padding:.58rem .68rem;border-bottom:1px solid #edf0ee;vertical-align:top;white-space:nowrap}tr:hover td{background:#faf7ef}details{margin:1rem 0}summary{cursor:pointer;font-weight:700}.file-links{display:flex;flex-wrap:wrap;gap:.5rem}.file-link{padding:.52rem .75rem;border:0;border-radius:999px;background:var(--teal-soft);font-size:.81rem;color:#145f5a}.file-link:hover,.file-link:focus-visible{background:#cbe6e1;outline:none;box-shadow:0 0 0 3px rgba(40,111,155,.2)}.unavailable{color:var(--muted);font-style:italic}a{color:#176c68}
footer{max-width:1540px;margin:auto;padding:1.5rem max(4vw,2rem);color:var(--muted);font-size:.75rem}
@media(max-width:950px){.mini-grid{grid-template-columns:1fr}.run-cards{grid-template-columns:1fr}.overview-strip{grid-template-columns:1fr 1fr}main{padding-inline:1rem}.plot-card{border-radius:17px}}
@media(max-width:560px){.overview-strip{grid-template-columns:1fr}header{padding-inline:1rem}nav{padding-inline:.7rem}.plot-card{padding:.75rem}.plot-actions{justify-content:flex-start}}
"""


def render_report(model: ReportModel, plots: list[dict[str, object]], output_path: Path) -> None:
    prepared_plots: list[dict[str, object]] = []
    for index, item in enumerate(plots, start=1):
        identity = "-".join(
            str(value)
            for value in (item.get("kind"), item.get("run"), item.get("modality"))
            if value not in {None, "", "ALL"}
        )
        readable_slug = re.sub(r"[^a-z0-9]+", "-", identity.lower()).strip("-") or "figure"
        identity_digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:8]
        slug = f"{readable_slug}-{identity_digest}"
        dom_id = f"qc-svg-{index}-{slug}"
        svg_download_id = f"qc-download-svg-{index}-{slug}"
        png_download_id = f"qc-download-png-{index}-{slug}"
        svg = str(item["svg"]).replace("<svg ", f'<svg id="{dom_id}" ', 1)
        prepared_plots.append(
            {
                **item,
                "svg": svg,
                "dom_id": dom_id,
                "svg_download_id": svg_download_id,
                "png_download_id": png_download_id,
                "svg_filename": f"{slug}.svg",
                "png_filename": f"{slug}.png",
            }
        )
    plots = prepared_plots
    retention = _retention_section(plots)
    mark = _single_plot_section(plots, "mark_composition", "DNA mark composition", "mark-composition")
    sample = _sample_section(plots)
    complexity = _dna_complexity_section(plots)
    rna_complexity = _rna_complexity_section(plots)
    document = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(model.library_name)} · TrESFlow QC</title><style>{CSS}</style></head>
<body><header><div class="brand">TrESFlow · sequencing QC</div><div class="hero"><h1>{esc(model.library_name)}</h1></div></header>
<nav><a href="#overview">Overview</a><a href="#retention">Retention</a>{'<a href="#mark-composition">DNA marks</a>' if mark else ''}{'<a href="#sample-composition">Sample barcodes</a>' if sample else ''}{'<a href="#complexity">DNA complexity</a>' if complexity else ''}{'<a href="#rna-complexity">RNA complexity</a>' if rna_complexity else ''}<a href="#details">Metrics</a></nav>
<main>{_overview(model)}{retention}{mark}{sample}{complexity}{rna_complexity}{_details(model)}</main>
<footer>Offline TrESFlow QC report · normalized schema 2.0 · no external resources</footer>
<script>(()=>{{
const prepareSvg=(svg)=>{{const clone=svg.cloneNode(true);clone.removeAttribute('id');clone.setAttribute('xmlns','http://www.w3.org/2000/svg');const viewBox=(clone.getAttribute('viewBox')||'0 0 1200 800').trim().split(/\\s+/).map(Number);clone.setAttribute('width',String(viewBox[2]));clone.setAttribute('height',String(viewBox[3]));return {{clone,viewBox,text:new XMLSerializer().serializeToString(clone)}};}};
const saveBlob=(blob,filename)=>{{const url=URL.createObjectURL(blob);const link=document.createElement('a');link.href=url;link.download=filename;document.body.appendChild(link);link.click();link.remove();setTimeout(()=>URL.revokeObjectURL(url),0);}};
document.querySelectorAll('[data-svg-download]').forEach(button=>{{button.addEventListener('click',()=>{{const svg=document.getElementById(button.dataset.svgDownload);if(!svg)return;const prepared=prepareSvg(svg);if(button.dataset.exportFormat==='svg'){{saveBlob(new Blob([prepared.text],{{type:'image/svg+xml;charset=utf-8'}}),button.dataset.filename);return;}}const source=new Blob([prepared.text],{{type:'image/svg+xml;charset=utf-8'}});const url=URL.createObjectURL(source);const image=new Image();image.onload=()=>{{const scale=2;const canvas=document.createElement('canvas');canvas.width=Math.ceil(prepared.viewBox[2]*scale);canvas.height=Math.ceil(prepared.viewBox[3]*scale);const context=canvas.getContext('2d');context.setTransform(scale,0,0,scale,0,0);context.fillStyle='#ffffff';context.fillRect(0,0,prepared.viewBox[2],prepared.viewBox[3]);context.drawImage(image,0,0,prepared.viewBox[2],prepared.viewBox[3]);URL.revokeObjectURL(url);canvas.toBlob(blob=>{{if(blob)saveBlob(blob,button.dataset.filename);}},'image/png');}};image.src=url;}});}});
}})();</script></body></html>'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8")
