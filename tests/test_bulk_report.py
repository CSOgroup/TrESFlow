"""Barcode abundance labels and denominators on the two bulk entry points."""
import sys
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from bulk_report import barcode_section, metric
from tresflow_qc.bulk_barcodes import BarcodePanel, abundance_svg, raw_panels, routed_panels


def test_raw_stage_denominators_and_groups():
    def row(kind, group, name, n, total):
        return dict(barcode_type=kind, group=group, category=name, count=n,
                    denominator_count=total, denominator_definition="measured pairs at this stage")
    observations = [row("sample_barcode", "a", "ATGC", 4, 5), row("sample_barcode", "__all__", "NoMatch", 1, 5),
                    row("dna_mark", "a", "target", 3, 4), row("dna_mark", "a", "NoMatch", 1, 4)]
    tag, split, mark = raw_panels("s", "dna", observations, {"ATGC": 8, "NoMatch": 2}, 10, 2)
    assert [p.denominator for p in (tag, split, mark)] == [10, 5, 4]
    assert tag.counts["ATGC", "a"] == 8 and split.counts["NoMatch", ""] == 1
    assert mark.scope == "s · DNA · a"
    assert "25.0%" in abundance_svg(mark)


def test_routed_marks_and_missing_tags():
    targets = [dict(sample="s", group="g", modality="dna", id=m, mark=m) for m in ("A", "B")]
    observations = {}
    for m, seq, n in [("A", "ATGC", 3), ("B", "CGTA", 7)]:
        observations["dna", m] = {**{tag: {"AAA": n} for tag in ("L1", "L2", "L3")},
                                  "SB": {"unavailable": n}, "MO": {seq: n}}
    panels = {p.tag: p for p in routed_panels(targets, observations)}
    assert panels["MO"].counts == {("ATGC", "A"): 3, ("CGTA", "B"): 7}
    assert all(p.denominator == 10 and p.stage == "Supplied routed pairs" for p in panels.values())
    svg = abundance_svg(panels["SB"])
    assert "Tag unavailable" in svg and "Unassigned" not in svg
    assert metric("Known count", 7, None, "raw pairs")[1:3] == (7, "unavailable")


def test_ranked_chart_keeps_unassigned_and_all_exact_counts():
    counts = Counter({(f"seq{i}", ""): 15-i for i in range(12)})
    counts["NoMatch", ""] = 1
    panel = BarcodePanel("s · RNA", "rna", "L1", "Split input", sum(counts.values()), "split-input pairs", counts)
    svg = abundance_svg(panel)
    ElementTree.fromstring(svg)
    assert "Top 8 of 12" in svg and "Unassigned" in svg and "seq11" not in svg
    html = barcode_section([panel])
    assert "seq11" in html.split("<details>")[1]
    assert html.count("data-svg-download=") == 2
