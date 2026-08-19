import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "lib"))

from tresflow_qc.cli import main as assessor_main
from tresflow_qc.core import (
    Branch,
    ReportModel,
    build_report_model,
    read_artifact_summary,
    read_duplicate_metrics,
    retention_display_stage,
)
from tresflow_qc.html_report import render_report
from tresflow_qc.plots import categorical_palette, generate_all_plots, generate_inline_plots


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def retention(path, modality, rows):
    lines = ["sample_id\tmodality\tgroup\tbranch\tmetric\tpairs\tunit"]
    for group, branch, metric, count in rows:
        lines.append(f"sample\t{modality}\t{group}\t{branch}\t{metric}\t{count}\tread_pairs")
    write(path, "\n".join(lines) + "\n")


def composition_row(modality, group, barcode_type, category, sequence, count, denominator):
    label = category
    definition = "ligation-accepted pairs" if barcode_type == "sample_barcode" else "sample-barcode-accepted pairs"
    return f"sample\t{modality}\t{group}\t{barcode_type}\t{category}\t{sequence}\t{label}\t{count}\t{count / denominator * 100:.6f}\t{denominator}\t{definition}\texisting decisions"


def synthetic_plot_model(
    *,
    run_count=1,
    modalities=("RNA", "DNA"),
    mark_count=2,
    barcode_counts=(25, 25, 25, 25),
    no_match_count=20,
    long_labels=False,
    roi=True,
    rna_saturation=True,
):
    branches = []
    retention_rows = []
    composition_rows = []
    complexity_rows = []
    qc_rows = []
    barcode_count = len(barcode_counts)
    sequences = [f"SB{index:03d}" for index in range(barcode_count)]
    marks = [
        (f"mark-{index + 1}" if not long_labels else f"mark-{index + 1}-with-a-deliberately-long-scientific-label")
        for index in range(mark_count)
    ]
    for run_index in range(run_count):
        run = f"run-{run_index + 1}" if not long_labels else f"independent-run-{run_index + 1}-with-a-deliberately-long-name"
        parent = f"parent-{run_index + 1}"
        group = f"group-{run_index + 1}"
        raw = 10_000 + run_index * 500
        if "RNA" in modalities:
            branches.append(Branch(run, parent, "RNA", "RNA", group, f"{parent}_{group}"))
            for order, (stage, count) in enumerate(
                (("Raw input", raw), ("After paired trimming", raw - 500), ("Called-cell final BAM", raw // 3)),
                start=1,
            ):
                retention_rows.append(
                    {
                        "run": run, "parent_sample": parent, "modality": "RNA", "branch": "RNA",
                        "stage_order": order, "stage": stage, "stage_scope": "shared",
                        "input_pairs": raw, "output_pairs": count, "retained_prev_pct": count / raw * 100,
                        "cumulative_raw_pct": count / raw * 100, "unit": "read pairs", "count_source": "synthetic", "subset_verified": "yes",
                    }
                )
            if rna_saturation:
                qc_rows.append(
                    {
                        "run": run, "parent_sample": parent, "modality": "RNA", "branch": "RNA",
                        "metric": "sequencing_saturation_pct", "value": 41.25 + run_index,
                        "unit": "percent",
                        "denominator": "STARsolo reads assigned to a unique GeneFull feature",
                        "source": f"{parent}.Solo.outGeneFull/Summary.csv",
                    }
                )
        if "DNA" in modalities:
            for mark_index, mark in enumerate(marks):
                branches.append(Branch(run, parent, "DNA", mark, group, f"{parent}_{group}_{mark_index}"))
                stages = (
                    ("Raw input", raw, "shared"),
                    ("After paired trimming", raw - 600, "shared"),
                    ("MO barcode accepted — all required barcodes accepted", raw - 2200, "shared"),
                    ("Mark-specific routing branch", max(1, raw // max(1, mark_count) - mark_index), "branch"),
                    ("Canonical NoDup final", max(1, raw // max(1, mark_count) - 100 - mark_index), "branch"),
                )
                for order, (stage, count, scope) in enumerate(stages, start=1):
                    retention_rows.append(
                        {
                            "run": run, "parent_sample": parent, "modality": "DNA", "branch": mark,
                            "stage_order": order, "stage": stage, "stage_scope": scope,
                            "input_pairs": raw, "output_pairs": count, "retained_prev_pct": count / raw * 100,
                            "cumulative_raw_pct": count / raw * 100, "unit": "read pairs", "count_source": "synthetic", "subset_verified": "yes",
                        }
                    )
                base_complexity = {
                    "run": run, "parent_sample": parent, "branch": mark,
                    "read_pairs_examined": 1000, "observed_unique_pairs": 940,
                    "pcr_or_library_duplicate_pairs": 35, "optical_duplicate_pairs": 25,
                    "observed_unique_fraction": 0.94, "pcr_or_library_duplicate_fraction": 0.035,
                    "optical_duplicate_fraction": 0.025, "percent_duplication": 6.0,
                    "estimated_library_size": 12_000, "source": "synthetic",
                }
                if roi:
                    for depth, yield_multiplier in ((1.0, 0.94), (2.0, 1.78), (4.0, 3.10), (8.0, 5.05), (16.0, 7.20)):
                        complexity_rows.append(
                            {
                                **base_complexity,
                                "roi_depth_multiplier": depth,
                                "roi_yield_multiplier": yield_multiplier,
                                "roi_estimated_unique_pairs": yield_multiplier * 1000,
                            }
                        )
                else:
                    complexity_rows.append(
                        {
                            **base_complexity,
                            "roi_depth_multiplier": None,
                            "roi_yield_multiplier": None,
                            "roi_estimated_unique_pairs": None,
                        }
                    )
        denominator = sum(barcode_counts) + no_match_count
        for modality in modalities:
            for sequence, count in zip(sequences, barcode_counts):
                label = sequence if not long_labels else f"barcode-label-{sequence}-with-a-deliberately-long-description"
                composition_rows.append(
                    {
                        "run": run, "parent_sample": parent, "modality": modality, "group": group,
                        "barcode_type": "sample_barcode", "category": sequence,
                        "barcode_sequence": sequence, "barcode_label": label, "count": count,
                        "percentage": count / denominator * 100 if denominator else 0.0,
                        "denominator_count": denominator, "denominator_definition": "synthetic eligible pairs",
                        "contract_level": "exact_sequence_counts_stats", "source": "synthetic",
                    }
                )
            composition_rows.append(
                {
                    "run": run, "parent_sample": parent, "modality": modality, "group": group,
                    "barcode_type": "sample_barcode", "category": "NoMatch", "barcode_sequence": "",
                    "barcode_label": "NoMatch", "count": no_match_count,
                    "percentage": no_match_count / denominator * 100 if denominator else 0.0,
                    "denominator_count": denominator, "denominator_definition": "synthetic eligible pairs",
                    "contract_level": "exact_sequence_counts_stats", "source": "synthetic:reads_without_bc",
                }
            )
        if "DNA" in modalities:
            mark_denominator = 10_000
            matched_total = mark_denominator - 23
            mark_values = [matched_total // mark_count] * mark_count
            mark_values[-1] += matched_total - sum(mark_values)
            for mark, count in zip(marks, mark_values):
                composition_rows.append(
                    {
                        "run": run, "parent_sample": parent, "modality": "DNA", "group": group,
                        "barcode_type": "dna_mark", "category": mark, "barcode_sequence": "",
                        "barcode_label": mark, "count": count, "percentage": count / mark_denominator * 100,
                        "denominator_count": mark_denominator, "denominator_definition": "synthetic mark-eligible pairs",
                        "contract_level": "exact_sequence", "source": "synthetic",
                    }
                )
            composition_rows.append(
                {
                    "run": run, "parent_sample": parent, "modality": "DNA", "group": group,
                    "barcode_type": "dna_mark", "category": "NoMatch", "barcode_sequence": "",
                    "barcode_label": "NoMatch", "count": 23, "percentage": 0.23,
                    "denominator_count": mark_denominator, "denominator_definition": "synthetic mark-eligible pairs",
                    "contract_level": "exact_sequence", "source": "synthetic",
                }
            )
    return ReportModel(
        library_name="synthetic", pipeline_version="test", generated_at="2026-01-01T00:00:00+00:00",
        branches=branches, retention=retention_rows, qc_metrics=qc_rows, barcode_composition=composition_rows,
        library_complexity=complexity_rows, insights=[], warnings=[], run_metadata={}, raw={},
    )


def expected_plot_set_from_fixture(model):
    """Independent test oracle derived directly from fixture rows."""
    expected = set()
    runs = sorted({branch.run for branch in model.branches})
    for run in runs:
        if any(row["run"] == run for row in model.retention):
            expected.add((run, "retention"))
        if any(
            row["run"] == run and row["barcode_type"] == "dna_mark"
            for row in model.barcode_composition
        ):
            expected.add((run, "mark_composition"))
        run_sources = {
            (branch.parent, branch.modality)
            for branch in model.branches
            if branch.run == run
        }
        if any(
            (row["parent_sample"], row["modality"]) in run_sources
            and row["barcode_type"] == "sample_barcode"
            and str(row.get("contract_level", "")).startswith("exact_sequence")
            for row in model.barcode_composition
        ):
            expected.add((run, "sample_composition"))
        if any(row["run"] == run for row in model.library_complexity):
            expected.add((run, "complexity"))
        if any(
            row["run"] == run and row["metric"] == "sequencing_saturation_pct"
            for row in model.qc_metrics
        ):
            expected.add((run, "rna_complexity"))
    return expected


def svg_geometry(svg):
    root = ET.fromstring(svg)
    view_box = [float(value) for value in root.attrib["viewBox"].split()]
    width, height = view_box[2], view_box[3]
    if float(root.attrib.get("width", 0)) != width or float(root.attrib.get("height", 0)) != height:
        raise AssertionError("SVG intrinsic dimensions must match its viewBox to preserve annotation size")
    layout_boxes = []
    categories = set()
    for element in root.iter():
        category = element.attrib.get("data-category")
        if category:
            categories.add(category)
        packed = element.attrib.get("data-layout-box")
        if packed:
            x, y, box_width, box_height = (float(value) for value in packed.split(","))
            layout_boxes.append((element.attrib.get("data-layout-scope", ""), x, y, box_width, box_height))
            if min(x, y, box_width, box_height) < -0.01 or x + box_width > width + 0.01 or y + box_height > height + 0.01:
                raise AssertionError(f"layout box outside viewBox: {(x, y, box_width, box_height)} vs {(width, height)}")
        tag = element.tag.rsplit("}", 1)[-1]
        if tag == "rect":
            x = float(element.attrib.get("x", 0)); y = float(element.attrib.get("y", 0))
            box_width = float(element.attrib.get("width", 0)); box_height = float(element.attrib.get("height", 0))
            if min(x, y, box_width, box_height) < -0.01 or x + box_width > width + 0.01 or y + box_height > height + 0.01:
                raise AssertionError(f"rectangle outside viewBox: {(x, y, box_width, box_height)} vs {(width, height)}")
        elif tag == "circle":
            x = float(element.attrib["cx"]); y = float(element.attrib["cy"]); radius = float(element.attrib.get("r", 0))
            if x - radius < -0.01 or y - radius < -0.01 or x + radius > width + 0.01 or y + radius > height + 0.01:
                raise AssertionError(f"circle outside viewBox: {(x, y, radius)} vs {(width, height)}")
        elif tag == "line":
            points = (
                (float(element.attrib.get("x1", 0)), float(element.attrib.get("y1", 0))),
                (float(element.attrib.get("x2", 0)), float(element.attrib.get("y2", 0))),
            )
            if any(x < -0.01 or y < -0.01 or x > width + 0.01 or y > height + 0.01 for x, y in points):
                raise AssertionError(f"line outside viewBox: {points} vs {(width, height)}")
        elif tag == "polyline":
            points = [tuple(float(value) for value in point.split(",")) for point in element.attrib.get("points", "").split()]
            if any(x < -0.01 or y < -0.01 or x > width + 0.01 or y > height + 0.01 for x, y in points):
                raise AssertionError(f"polyline outside viewBox: {points} vs {(width, height)}")
        elif tag == "path":
            current_x = current_y = 0.0
            tokens = re.findall(r"[MHV]|-?[0-9]+(?:\.[0-9]+)?", element.attrib.get("d", ""))
            index = 0
            while index < len(tokens):
                command = tokens[index]; index += 1
                if command == "M":
                    current_x = float(tokens[index]); current_y = float(tokens[index + 1]); index += 2
                elif command == "H":
                    current_x = float(tokens[index]); index += 1
                elif command == "V":
                    current_y = float(tokens[index]); index += 1
                if current_x < -0.01 or current_y < -0.01 or current_x > width + 0.01 or current_y > height + 0.01:
                    raise AssertionError(f"path outside viewBox: {(current_x, current_y)} vs {(width, height)}")
    return width, height, layout_boxes, categories


def assert_non_overlapping_layout(test_case, boxes):
    by_scope = {}
    for scope, x, y, width, height in boxes:
        by_scope.setdefault(scope, []).append((x, y, width, height))
    for scope, scoped_boxes in by_scope.items():
        for left_index, left in enumerate(scoped_boxes):
            for right in scoped_boxes[left_index + 1:]:
                overlap_width = min(left[0] + left[2], right[0] + right[2]) - max(left[0], right[0])
                overlap_height = min(left[1] + left[3], right[1] + right[3]) - max(left[1], right[1])
                test_case.assertFalse(
                    overlap_width > 0.5 and overlap_height > 0.5,
                    f"overlap in {scope}: {left} vs {right}",
                )


class QcFixture:
    def __init__(self, root, *, rna=True, dna=True, dt=True, roi=True, exact=True, groups=None, marks=None, zero_no_match=False):
        self.root = root
        stats = root / "TrES_Stats"
        contract = root / "pipeline_info" / "derived_contract"
        groups = groups or {"g1": ["AAA"], "g2": ["CCC"]}
        marks = marks or {"g1": ["H3A", "H3B"], "g2": ["H3C"]}
        configured_barcodes = [barcode for values in groups.values() for barcode in values]

        def tag_counts(total_assigned, no_match):
            base = total_assigned // len(configured_barcodes)
            values = [base] * len(configured_barcodes)
            values[-1] += total_assigned - sum(values)
            return "\n".join(
                [f"{count}\t{barcode}" for barcode, count in zip(configured_barcodes, values)]
                + [f"{no_match}\tNoMatch"]
            )
        sb_lines = ["sample\tsb_group\tsb_bc"]
        for group, barcodes in groups.items():
            for barcode in barcodes:
                sb_lines.append(f"sample\t{group}\t{barcode}")
        if rna:
            write(contract / "rna_sb_group_map.tsv", "\n".join(sb_lines) + "\n")
        if dna:
            write(contract / "dna_sb_group_map.tsv", "\n".join(sb_lines) + "\n")
            mo_lines = ["sample\tsb_group\tmark\tmo_bc"]
            for group, group_marks in marks.items():
                for index, mark in enumerate(group_marks):
                    mo_lines.append(f"sample\t{group}\t{mark}\tMO{index}{group}")
            write(contract / "dna_mo_map.tsv", "\n".join(mo_lines) + "\n")
        if rna:
            rna_tagged = 100 if zero_no_match else 40
            rna_no_match = 0 if zero_no_match else 60
            write(stats / "sample.rna_sample_barcode.counts.tsv", tag_counts(rna_tagged, rna_no_match))
            write(stats / "sample.rna_sample_barcode.stats.tsv", f"reads\t100\t100\nbc_reads\t{rna_tagged}\t{rna_tagged}\nreads_without_bc\t{rna_no_match}\t{rna_no_match}\n")
            for layer in (1, 2, 3):
                write(stats / f"sample.rna_cell.stats_L{layer}.tsv", f"reads_with_L{layer}\t70\t70\n")
            rna_split_rows = [("__all__", "__all__", "split_input_pairs", 90), ("__all__", "__all__", "joint_barcode_accepted_pairs", 40)]
            for group in groups:
                rna_split_rows += [(group, "__all__", "routed_group_pairs", 20), (group, "RNA", "routed_branch_pairs", 20)]
                split = f"sample_{group}"
                retention(stats / f"{split}.rna_filter_retention.tsv", "rna", [(group, "RNA", "paired_filter_pairs", 18), (group, "RNA", "canonical_pairs", 16), (group, "RNA", "called_cell_pairs", 14)])
                write(
                    root / "rna_align" / f"{split}.Solo.outGeneFull" / "Summary.csv",
                    "Number of Reads,20\nSequencing Saturation,0.375\nEstimated Number of Cells,2\nMedian UMI per Cell,4\n",
                )
            retention(stats / "sample.rna_read_retention.tsv", "rna", rna_split_rows)
            gate_header = "sample_id\tmodality\tmetric\tpairs\tunit\tdenominator_metric\tdenominator_pairs\tdefinition\tsource\n"
            if exact:
                rna_accepted = 60 if zero_no_match else 40
                write(stats / "sample.rna_barcode_gates.tsv", gate_header + "sample\trna\tsplit_input_pairs\t90\tread_pairs\tsplit_input_pairs\t90\td\ts\n" + "sample\trna\tligation_barcode_accepted_pairs\t60\tread_pairs\tsplit_input_pairs\t90\td\ts\n" + f"sample\trna\tsample_barcode_accepted_pairs\t{rna_accepted}\tread_pairs\tsplit_input_pairs\t90\td\ts\n")
        if dna:
            dna_tagged = 100 if zero_no_match else 50
            dna_no_match = 0 if zero_no_match else 50
            write(stats / "sample.dna_sample_barcode.counts.tsv", tag_counts(dna_tagged, dna_no_match))
            write(stats / "sample.dna_sample_barcode.stats.tsv", f"reads\t100\t100\nbc_reads\t{dna_tagged}\t{dna_tagged}\nreads_without_bc\t{dna_no_match}\t{dna_no_match}\n")
            for layer in (1, 2, 3):
                write(stats / f"sample.dna_cell.stats_L{layer}.tsv", f"reads_with_L{layer}\t65\t65\n")
            write(stats / "sample.dna_modality.stats.tsv", "reads\t100\t100\nbc_reads\t40\t40\n")
            dna_split_rows = [("__all__", "__all__", "split_input_pairs", 70), ("__all__", "__all__", "joint_barcode_accepted_pairs", 40)]
            for group, group_marks in marks.items():
                for mark in group_marks:
                    dna_split_rows.append((group, mark, "routed_branch_pairs", 10))
                    split = f"sample_{group}_{mark}"
                    retention(stats / f"{split}.dna_alignment_retention.tsv", "dna", [(group, mark, "proper_pair_primary_pairs", 9)])
                    write(stats / "qc" / "samtools" / f"dna.{split}.markeddup.flagstat", "8 + 0 read1\n")
                    write(stats / "qc" / "samtools" / f"dna.{split}.nodup.flagstat", "7 + 0 read1\n")
                    histogram = "\n## HISTOGRAM\tjava.lang.Double\nBIN\tCoverageMult\tall_sets\n1.0\t1.05\t8\n2.0\t1.90\t1\n5.0\t4.10\t0\n10.0\t7.20\t0\n" if roi else ""
                    write(root / "dna_align" / f"{split}.DuplicateMetrics.txt", "## METRICS CLASS\tpicard.sam.DuplicationMetrics\nLIBRARY\tREAD_PAIRS_EXAMINED\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\nlib\t9\t2\t1\t0.222222\t30\n" + histogram)
            retention(stats / "sample.dna_read_retention.tsv", "dna", dna_split_rows)
            if exact:
                gate_header = "sample_id\tmodality\tmetric\tpairs\tunit\tdenominator_metric\tdenominator_pairs\tdefinition\tsource\n"
                dna_accepted = 55 if zero_no_match else 50
                write(stats / "sample.dna_barcode_gates.tsv", gate_header + "sample\tdna\tsplit_input_pairs\t70\tread_pairs\tsplit_input_pairs\t70\td\ts\n" + "sample\tdna\tligation_barcode_accepted_pairs\t55\tread_pairs\tsplit_input_pairs\t70\td\ts\n" + f"sample\tdna\tsample_barcode_accepted_pairs\t{dna_accepted}\tread_pairs\tsplit_input_pairs\t70\td\ts\n" + f"sample\tdna\tmodality_barcode_accepted_pairs\t{dna_accepted}\tread_pairs\tsplit_input_pairs\t70\td\ts\n")
            if dt:
                write(stats / "sample.dual_tag_artifact_filter.summary.tsv", "sample_id\ttagmentation\tinput_pairs\tretained_pairs\trejected_pairs\trejected_fraction\tr1_with_signature\tr2_with_signature\tcutadapt_version\tsignature_fasta_sha256\n" + "sample\tdual\t80\t70\t10\t0.125\t8\t5\t5.0\t67c6f1789ef5e36492562203ac38fc13fa901058047ed2bd37b304d85a30ae0f\n")
        if exact:
            header = "sample_id\tmodality\tgroup\tbarcode_type\tcategory\tbarcode_sequence\tbarcode_label\tcount\tpercentage\tdenominator_count\tdenominator_definition\tsource"
            barcodes = [(group, barcode) for group, values in groups.items() for barcode in values]
            if rna:
                accepted = 40
                base = accepted // len(barcodes)
                counts = [base] * len(barcodes)
                counts[-1] += accepted - sum(counts)
                no_match = 0 if zero_no_match else 20
                if zero_no_match:
                    counts[-1] += 20
                rows = [composition_row("rna", group, "sample_barcode", barcode, barcode, count, 60) for (group, barcode), count in zip(barcodes, counts)]
                rows.append(composition_row("rna", "__all__", "sample_barcode", "NoMatch", "", no_match, 60))
                write(stats / "sample.rna_barcode_composition.tsv", header + "\n" + "\n".join(rows) + "\n")
            if dna:
                sample_counts = [30, 20] if len(barcodes) == 2 else [50 // len(barcodes)] * len(barcodes)
                sample_counts[-1] += 50 - sum(sample_counts)
                sample_no_match = 0 if zero_no_match else 5
                if zero_no_match:
                    sample_counts[-1] += 5
                rows = [composition_row("dna", group, "sample_barcode", barcode, barcode, count, 55) for (group, barcode), count in zip(barcodes, sample_counts)]
                rows.append(composition_row("dna", "__all__", "sample_barcode", "NoMatch", "", sample_no_match, 55))
                for group, group_marks in marks.items():
                    denominator = sum(count for (barcode_group, _), count in zip(barcodes, sample_counts) if barcode_group == group)
                    mark_total = denominator
                    mark_values = [mark_total // len(group_marks)] * len(group_marks)
                    mark_values[-1] += mark_total - sum(mark_values)
                    rows += [composition_row("dna", group, "dna_mark", mark, "", count, denominator) for mark, count in zip(group_marks, mark_values)]
                    rows.append(composition_row("dna", group, "dna_mark", "NoMatch", "", denominator - mark_total, denominator))
                write(stats / "sample.dna_barcode_composition.tsv", header + "\n" + "\n".join(rows) + "\n")


class SharedQcTests(unittest.TestCase):
    def test_modern_contract_stages_composition_dt_and_complexity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root)
            model = build_report_model(root, "fixture", "1.2.3")
            rna = [row["stage"] for row in model.retention if row["run"] == "g1" and row["modality"] == "RNA"]
            self.assertEqual(rna, ["Raw input", "After paired trimming", "L1–L2–L3 accepted", "Sample barcode accepted — all required barcodes accepted", "Properly paired", "Canonical chromosomes", "Called-cell final BAM"])
            dna = [row["stage"] for row in model.retention if row["run"] == "g1" and row["branch"] == "H3A"]
            self.assertEqual(dna, ["Raw input", "After paired trimming", "After DT artifact filter", "L1–L2–L3 accepted", "Sample barcode accepted", "MO barcode accepted — all required barcodes accepted", "Mark-specific routing branch", "Properly paired after blacklist", "Canonical marked-duplicate BAM", "Canonical NoDup final"])
            removed = {"Routed to run", "Entered STAR", "STAR mapped primary", "BWA primary", "After blacklist removal", "Mapped after blacklist", "Joint barcode accepted"}
            self.assertFalse(removed & {row["stage"] for row in model.retention})
            sample_rows = [row for row in model.barcode_composition if row["barcode_type"] == "sample_barcode" and row["modality"] == "DNA"]
            self.assertEqual({row["barcode_sequence"] for row in sample_rows}, {"AAA", "CCC", ""})
            self.assertEqual(sum(row["count"] for row in sample_rows), 100)
            self.assertTrue(all(row["contract_level"] == "exact_sequence_counts_stats" for row in sample_rows))
            no_match = next(row for row in sample_rows if row["category"] == "NoMatch")
            self.assertEqual(no_match["count"], 50)
            self.assertIn("reads_without_bc", no_match["source"])
            dt = {row["metric"]: row["value"] for row in model.qc_metrics if str(row["metric"]).startswith("dual_tag_filter_")}
            self.assertEqual((dt["dual_tag_filter_both_mates_signature_pairs"], dt["dual_tag_filter_r1_only_signature_pairs"], dt["dual_tag_filter_r2_only_signature_pairs"]), (3, 5, 2))
            row = next(row for row in model.library_complexity if row["run"] == "g1" and row["branch"] == "H3A")
            self.assertEqual(row["observed_unique_pairs"] + row["pcr_or_library_duplicate_pairs"] + row["optical_duplicate_pairs"], row["read_pairs_examined"])

    def test_dt_reconciliation_errors_and_non_dt_omits_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=False, dt=False, groups={"solo": ["AAA"]}, marks={"solo": ["M"]})
            model = build_report_model(root)
            self.assertNotIn("After DT artifact filter", {row["stage"] for row in model.retention})
            self.assertFalse(any(str(row["metric"]).startswith("dual_tag_filter_") for row in model.qc_metrics))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.summary.tsv"
            write(path, "sample_id\tinput_pairs\tretained_pairs\trejected_pairs\tr1_with_signature\tr2_with_signature\nS\t10\t8\t3\t2\t2\n")
            with self.assertRaisesRegex(ValueError, "accounting mismatch"):
                read_artifact_summary(path)

    def test_legacy_uses_coarse_stage_without_fabrication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, exact=False)
            model = build_report_model(root)
            stages = {row["stage"] for row in model.retention}
            self.assertIn("Joint barcode accepted", stages)
            self.assertNotIn("L1–L2–L3 accepted", stages)
            self.assertTrue(any(row["barcode_type"] == "sample_barcode" for row in model.barcode_composition))
            self.assertTrue(any("coarse" in warning for warning in model.warnings))

    def test_grouped_legacy_sample_composition_is_retained_in_table_but_not_plotted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "input"
            QcFixture(root, rna=True, dna=False, groups={"solo": ["AAA"]})
            path = root / "TrES_Stats" / "sample.rna_barcode_composition.tsv"
            (root / "TrES_Stats" / "sample.rna_sample_barcode.counts.tsv").unlink()
            write(
                path,
                "sample_id\tmodality\tgroup\tbarcode_type\tcategory\tcount\tpercentage\tdenominator_count\tdenominator_definition\tsource\n"
                "sample\trna\t__all__\tsample_barcode\tsolo\t40\t66.666667\t60\tligation-accepted pairs\tlegacy\n"
                "sample\trna\t__all__\tsample_barcode\tNoMatch\t20\t33.333333\t60\tligation-accepted pairs\tlegacy\n",
            )
            model = build_report_model(root)
            self.assertEqual({row["contract_level"] for row in model.barcode_composition}, {"legacy_group"})
            plots = generate_all_plots(model, Path(tmp) / "plots")
            self.assertFalse(any(item["kind"] == "sample_composition" for item in plots))
            self.assertTrue(any("per-barcode composition unavailable" in warning for warning in model.warnings))

    def test_picard_roi_and_missing_roi(self):
        for roi, expected_points in ((True, 4), (False, 1)):
            with self.subTest(roi=roi), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                QcFixture(root, rna=False, roi=roi, groups={"solo": ["AAA"]}, marks={"solo": ["M"]})
                model = build_report_model(root)
                rows = [row for row in model.library_complexity if row["branch"] == "M"]
                self.assertEqual(len(rows), expected_points)
                self.assertEqual(rows[0]["optical_duplicate_pairs"], 1)
                if roi:
                    self.assertEqual([row["roi_depth_multiplier"] for row in rows], [1.0, 2.0, 5.0, 10.0])
                else:
                    self.assertIsNone(rows[0]["roi_depth_multiplier"])

    def test_zero_examined_picard_metrics_render_without_division_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=False, groups={"solo": ["AAA"]}, marks={"solo": ["M"]}, roi=False)
            metrics = root / "dna_align" / "sample_solo_M.DuplicateMetrics.txt"
            write(
                metrics,
                "## METRICS CLASS\tpicard.sam.DuplicationMetrics\n"
                "LIBRARY\tREAD_PAIRS_EXAMINED\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n"
                "lib\t0\t0\t0\t0\t0\n",
            )
            model = build_report_model(root)
            plots = generate_all_plots(model, root / "plots")
            self.assertTrue(any(item["kind"] == "complexity" for item in plots))

    def test_zero_nomatch_and_arbitrary_long_categories(self):
        groups = {"a very long group name": [f"BC{i:02d}" for i in range(14)]}
        marks = {"a very long group name": ["mark alpha with long name", "mark beta"]}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=False, groups=groups, marks=marks, zero_no_match=True)
            model = build_report_model(root)
            no_match = [row for row in model.barcode_composition if row["category"] == "NoMatch"]
            self.assertTrue(no_match)
            self.assertTrue(all(row["count"] == 0 for row in no_match if row["barcode_type"] == "sample_barcode"))
            self.assertEqual(len([row for row in model.barcode_composition if row["barcode_type"] == "sample_barcode"]), 15)

    def test_rna_only_and_dna_only(self):
        for rna, dna in ((True, False), (False, True)):
            with self.subTest(rna=rna, dna=dna), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                QcFixture(root, rna=rna, dna=dna, groups={"solo": ["AAA"]}, marks={"solo": ["M"]})
                model = build_report_model(root)
                self.assertEqual({branch.modality for branch in model.branches}, {"RNA" if rna else "DNA"})

    def test_parent_wide_sample_source_is_attached_only_to_applicable_independent_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(
                root,
                rna=True,
                dna=False,
                groups={"run-a": ["BAR-A"], "run-b": ["BAR-B"]},
            )
            model = build_report_model(root)
            plots = generate_inline_plots(model)
            sample_plots = [item for item in plots if item["kind"] == "sample_composition"]
            applicable = {branch.run for branch in model.branches if branch.modality == "RNA"}
            self.assertEqual({item["run"] for item in sample_plots}, applicable)
            source_categories = {
                str(row["barcode_sequence"] or row["category"])
                for row in model.barcode_composition
                if row["barcode_type"] == "sample_barcode"
            }
            for plot in sample_plots:
                _, _, _, represented = svg_geometry(plot["svg"])
                self.assertEqual(represented, source_categories)

    def test_cli_tables_figures_and_offline_report_dom(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "input"
            output = Path(tmp) / "output"
            QcFixture(root, groups={"solo": ["AAA"]}, marks={"solo": ["M"]})
            old_mpl = os.environ.get("MPLCONFIGDIR")
            os.environ["MPLCONFIGDIR"] = str(Path(tmp) / "mpl")
            try:
                self.assertEqual(assessor_main([str(root), "--output-dir", str(output), "--library-name", "QC fixture"]), 0)
            finally:
                if old_mpl is None:
                    os.environ.pop("MPLCONFIGDIR", None)
                else:
                    os.environ["MPLCONFIGDIR"] = old_mpl
            for name in ("read_retention.tsv", "qc_metrics.tsv", "barcode_composition.tsv", "library_complexity.tsv", "tres_report_metrics.json", "tres_report.html", "solo_read_retention.svg"):
                self.assertTrue((output / name).is_file(), name)
            document = (output / "tres_report.html").read_text(encoding="utf-8")
            for section in ("overview", "retention", "mark-composition", "sample-composition", "complexity", "details"):
                self.assertIn(f'id="{section}"', document)
            for removed in ("Automated QC summary", "Factual observations", "Exact gates and assigned categories", "Methods and provenance"):
                self.assertNotIn(removed, document)
            self.assertIn("<svg", document)
            self.assertNotIn("data:image/png;base64,", document)
            self.assertNotIn("cdn.", document.lower())
            payload = json.loads((output / "tres_report_metrics.json").read_text())
            self.assertEqual(payload["schema_version"], "2.0")

    def test_pipeline_contract_is_html_plus_four_tables_with_inline_svg(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "arbitrary-pipeline-output"
            output = Path(tmp) / "output"
            QcFixture(root, groups={"solo": ["AAA"]}, marks={"solo": ["M"]})
            code = assessor_main(
                [
                    str(root), "--output-dir", str(output), "--pipeline-version", "v8.4.2+3.g1234567",
                    "--no-json", "--no-standalone-figures",
                ]
            )
            self.assertEqual(code, 0)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {"tres_report.html", "read_retention.tsv", "qc_metrics.tsv", "barcode_composition.tsv", "library_complexity.tsv"},
            )
            document = (output / "tres_report.html").read_text(encoding="utf-8")
            self.assertIn("<title>arbitrary-pipeline-output · TrESFlow QC</title>", document)
            self.assertIn("<h1>arbitrary-pipeline-output</h1>", document)
            self.assertNotIn("<h1>output</h1>", document)
            self.assertIn("<strong>v8.4.2+3.g1234567</strong>", document)
            self.assertIn("<svg", document)
            self.assertNotIn("Exact read-pair accounting, barcode composition and DNA library complexity", document)
            self.assertNotIn("Independent run / group", document)
            self.assertNotIn('<div class="eyebrow">Independent run</div>', document)
            self.assertIn("<h3>solo</h3>", document)
            self.assertIn("100 / 14 (14.0% of raw)", document)
            self.assertIn("100 / 7 (7.0% of raw)", document)
            self.assertIn("RNA fully barcode-assigned pairs", document)
            self.assertIn("40 (40.0% of raw)", document)
            self.assertIn("DNA fully barcode-assigned pairs", document)
            self.assertIn("50 (50.0% of raw)", document)
            self.assertIn("DNA duplicate components", document)
            self.assertIn("PCR/library:", document)
            self.assertIn("Optical:", document)
            self.assertIn("RNA sequencing saturation", document)
            self.assertNotIn("Observed DNA duplicates", document)
            self.assertNotIn("RNA duplicates", document)
            mark_section = document.split('id="mark-composition"', 1)[1].split('id="sample-composition"', 1)[0]
            self.assertNotIn("NoMatch", mark_section)
            self.assertIn("NoMatch", document.split('id="sample-composition"', 1)[1])
            ids = re.findall(r'\bid="([^"]+)"', document)
            self.assertEqual(len(ids), len(set(ids)))
            self.assertFalse(any(link.startswith(("http://", "https://", "//")) for link in re.findall(r'href="([^"]+)"', document)))
            svg_ids = set(re.findall(r'<svg\s+id="([^"]+)"', document))
            download_targets = set(re.findall(r'data-svg-download="([^"]+)"', document))
            self.assertEqual(svg_ids, download_targets)
            self.assertEqual(document.count("Download SVG"), len(svg_ids))
            self.assertEqual(document.count("Download PNG"), len(svg_ids))

    def test_report_title_default_tracks_assessed_directory_and_override_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "results-alpha"
            second = Path(tmp) / "results-beta"
            QcFixture(first, rna=True, dna=False, groups={"cohort-a": ["AAA"]})
            QcFixture(second, rna=True, dna=False, groups={"cohort-b": ["CCC"]})
            for root in (first, second):
                output = Path(tmp) / f"render-{root.name}"
                self.assertEqual(
                    assessor_main([str(root), "--output-dir", str(output), "--no-json", "--no-standalone-figures"]),
                    0,
                )
                document = (output / "tres_report.html").read_text(encoding="utf-8")
                self.assertIn(f"<title>{root.name} · TrESFlow QC</title>", document)
                self.assertIn(f"<h1>{root.name}</h1>", document)
                self.assertNotIn(f"<h1>{output.name}</h1>", document)

            override_output = Path(tmp) / "override-render"
            self.assertEqual(
                assessor_main(
                    [
                        str(first), "--output-dir", str(override_output), "--title", "Explicit synthetic title",
                        "--no-json", "--no-standalone-figures",
                    ]
                ),
                0,
            )
            override_document = (override_output / "tres_report.html").read_text(encoding="utf-8")
            self.assertIn("<h1>Explicit synthetic title</h1>", override_document)

    def test_release_version_resolver_uses_exact_tag_development_baseline_and_manifest_fallback(self):
        resolver = REPO_ROOT / "bin" / "resolve_tresflow_release_version.sh"
        manifest_text = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")
        manifest_match = re.search(r"(?m)^\s*version\s*=\s*'([^']+)'", manifest_text)
        self.assertIsNotNone(manifest_match)
        self.assertRegex(manifest_match.group(1), r"^v[0-9]+\.[0-9]+\.[0-9]+(?:[.-][0-9A-Za-z]+)*$")
        with tempfile.TemporaryDirectory() as tmp:
            repository = Path(tmp) / "synthetic-repository"
            repository.mkdir()
            subprocess.run(["git", "init", "-q", str(repository)], check=True)
            subprocess.run(["git", "-C", str(repository), "config", "user.email", "fixture@example.invalid"], check=True)
            subprocess.run(["git", "-C", str(repository), "config", "user.name", "Fixture"], check=True)
            write(repository / "tracked.txt", "first\n")
            subprocess.run(["git", "-C", str(repository), "add", "tracked.txt"], check=True)
            subprocess.run(["git", "-C", str(repository), "commit", "-q", "-m", "first"], check=True)
            release_tag = "v8.4.2"
            subprocess.run(["git", "-C", str(repository), "tag", release_tag], check=True)
            exact = subprocess.run(
                ["bash", str(resolver), str(repository), "3.2.1dev"],
                check=True, capture_output=True, text=True,
            ).stdout.strip()
            self.assertEqual(exact, release_tag)

            write(repository / "tracked.txt", "second\n")
            subprocess.run(["git", "-C", str(repository), "add", "tracked.txt"], check=True)
            subprocess.run(["git", "-C", str(repository), "commit", "-q", "-m", "second"], check=True)
            development = subprocess.run(
                ["bash", str(resolver), str(repository), "3.2.1dev"],
                check=True, capture_output=True, text=True,
            ).stdout.strip()
            self.assertRegex(development, rf"^{re.escape(release_tag)}\+1\.g[0-9a-f]{{7}}$")

            archive = Path(tmp) / "source-archive"
            archive.mkdir()
            fallback = subprocess.run(
                ["bash", str(resolver), str(archive), "3.2.1dev"],
                check=True, capture_output=True, text=True,
            ).stdout.strip()
            self.assertEqual(fallback, "v3.2.1dev")

    def test_overview_card_titles_use_explicit_groups_for_reordered_long_unicode_metadata(self):
        model = synthetic_plot_model(
            run_count=2,
            modalities=("RNA", "DNA"),
            mark_count=2,
            barcode_counts=(5, 3, 2, 1),
        )
        group_for_run = {
            "run-1": "cohort-with-a-deliberately-long-samplesheet-group-name",
            "run-2": "ομάδα-測試",
        }
        model.branches = [
            Branch(
                branch.run,
                branch.parent,
                branch.modality,
                branch.branch,
                group_for_run[branch.run],
                branch.split_id,
            )
            for branch in reversed(model.branches)
        ]
        model.run_metadata = {
            "samples": [
                {"id": branch.parent, "modality": branch.modality.lower(), "groups": [branch.group]}
                for branch in reversed(model.branches)
            ]
        }
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report.html"
            render_report(model, generate_inline_plots(model), report)
            document = report.read_text(encoding="utf-8")
        for group in group_for_run.values():
            self.assertIn(f"<h3>{group}</h3>", document)
        self.assertNotIn('<div class="eyebrow">Independent run</div>', document)
        self.assertNotIn("<h3>run-1</h3>", document.split('id="retention"', 1)[0])
        self.assertNotIn("<h3>run-2</h3>", document.split('id="retention"', 1)[0])

    def test_retention_mapped_label_is_display_only_and_source_remains_explicit(self):
        model = synthetic_plot_model(
            run_count=1,
            modalities=("RNA", "DNA"),
            mark_count=1,
            barcode_counts=(7,),
        )
        rna = next(row for row in model.retention if row["modality"] == "RNA" and row["stage"] == "Called-cell final BAM")
        dna = next(row for row in model.retention if row["modality"] == "DNA" and row["stage"] == "Canonical NoDup final")
        rna["stage"] = "Properly paired"
        dna["stage"] = "Properly paired after blacklist"
        original = [(row["stage"], row["output_pairs"], row["count_source"]) for row in model.retention]
        retention_plot = next(item for item in generate_inline_plots(model) if item["kind"] == "retention")
        svg = retention_plot["svg"]
        root = ET.fromstring(svg)
        stage_text = [
            "".join(element.itertext())
            for element in root.iter()
            if element.tag.rsplit("}", 1)[-1] == "text" and element.attrib.get("class") == "stage"
        ]
        self.assertEqual(stage_text.count("MappedMapped"), 2)
        tooltips = [element.text or "" for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "title"]
        self.assertTrue(any(text.startswith("Mapped:") and 'underlying accounting stage "Properly paired"' in text for text in tooltips))
        self.assertTrue(any(text.startswith("Mapped:") and 'underlying accounting stage "Properly paired after blacklist"' in text for text in tooltips))
        self.assertEqual(
            [(row["stage"], row["output_pairs"], row["count_source"]) for row in model.retention],
            original,
        )
        self.assertEqual(retention_display_stage("Properly paired"), "Mapped")
        self.assertEqual(retention_display_stage("Properly paired after blacklist"), "Mapped")

        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "report.html"
            render_report(model, generate_inline_plots(model), report)
            document = report.read_text(encoding="utf-8")
        table = document.split("<summary>Read-pair retention</summary>", 1)[1].split("</details>", 1)[0]
        self.assertEqual(table.count("<td>Mapped</td>"), 2)
        self.assertIn("<td>Properly paired</td>", table)
        self.assertIn("<td>Properly paired after blacklist</td>", table)

    def test_sample_composition_has_external_key_only_for_every_category(self):
        for counts in ((9,), (9, 2, 1, 0), tuple([1] * 18 + [0] * 12)):
            with self.subTest(category_count=len(counts)):
                model = synthetic_plot_model(
                    run_count=1,
                    modalities=("RNA", "DNA"),
                    mark_count=2,
                    barcode_counts=counts,
                    no_match_count=300,
                    long_labels=True,
                )
                plots = generate_inline_plots(model)
                sample = next(item for item in plots if item["kind"] == "sample_composition")
                root = ET.fromstring(sample["svg"])
                self.assertFalse(any("data-in-bar-label" in element.attrib for element in root.iter()))
                bars = [element for element in root.iter() if element.attrib.get("data-segment-bar") == "true"]
                self.assertTrue(bars)
                self.assertTrue(all(not any(child.tag.rsplit("}", 1)[-1] == "text" for child in element) for element in bars))
                external = {
                    element.attrib["data-static-category"]
                    for element in root.iter()
                    if "data-static-category" in element.attrib
                }
                expected = {
                    str(row["barcode_sequence"] or row["category"])
                    for row in model.barcode_composition
                    if row["barcode_type"] == "sample_barcode"
                }
                self.assertEqual(external, expected)
                _, _, boxes, _ = svg_geometry(sample["svg"])
                assert_non_overlapping_layout(self, boxes)

                mark = next(item for item in plots if item["kind"] == "mark_composition")
                self.assertTrue(any("data-in-bar-label" in element.attrib for element in ET.fromstring(mark["svg"]).iter()))

    def test_sample_barcode_counts_and_stats_must_reconcile(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=True, dna=False, groups={"solo": ["AAA"]})
            write(root / "TrES_Stats" / "sample.rna_sample_barcode.stats.tsv", "reads\t100\t100\nbc_reads\t40\t40\nreads_without_bc\t59\t59\n")
            with self.assertRaisesRegex(ValueError, "NoMatch disagreement"):
                build_report_model(root)

    def test_configured_barcode_absent_from_counts_is_preserved_as_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=True, dna=False, groups={"run-a": ["BAR-A", "BAR-B"]})
            write(
                root / "TrES_Stats" / "sample.rna_sample_barcode.counts.tsv",
                "40\tBAR-A\n60\tNoMatch\n",
            )
            model = build_report_model(root)
            sample_rows = [
                row for row in model.barcode_composition
                if row["barcode_type"] == "sample_barcode" and row["modality"] == "RNA"
            ]
            independently_expected = {"BAR-A": 40, "BAR-B": 0, "NoMatch": 60}
            self.assertEqual(
                {str(row["barcode_sequence"] or row["category"]): int(row["count"]) for row in sample_rows},
                independently_expected,
            )
            self.assertEqual(sum(independently_expected.values()), int(sample_rows[0]["denominator_count"]))

    def test_missing_optional_inputs_do_not_create_empty_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "input"
            QcFixture(root, rna=True, dna=False, groups={"solo": ["AAA"]})
            model = build_report_model(root)
            output = Path(tmp) / "plots"
            with mock.patch("tresflow_qc.plots.generate_complexity_plots", wraps=generate_all_plots.__globals__["generate_complexity_plots"]):
                plots = generate_all_plots(model, output)
            self.assertFalse(any(item["kind"] in {"complexity", "mark_composition"} for item in plots))

    def test_adaptive_plot_matrix_is_deterministic_and_within_viewboxes(self):
        scenarios = [
            {
                "run_count": 1, "modalities": ("RNA",), "mark_count": 1,
                "barcode_counts": (0,), "no_match_count": 10_000, "roi": False,
            },
            {
                "run_count": 2, "modalities": ("DNA",), "mark_count": 2,
                "barcode_counts": (100, 100, 100, 100), "no_match_count": 0, "roi": True,
            },
            {
                "run_count": 4, "modalities": ("RNA", "DNA"), "mark_count": 6,
                "barcode_counts": tuple([1] * 7 + [140] * 5), "no_match_count": 2_500,
                "long_labels": True, "roi": False,
            },
            {
                "run_count": 2, "modalities": ("RNA", "DNA"), "mark_count": 1,
                "barcode_counts": tuple([0, 1, 1, 1, 1] + [75] * 25), "no_match_count": 10,
                "long_labels": True, "roi": True,
            },
        ]
        observed_run_counts = set()
        observed_mark_counts = set()
        observed_barcode_counts = set()
        sample_widths = {}
        retention_heights = {}
        for scenario in scenarios:
            with self.subTest(scenario=scenario):
                model = synthetic_plot_model(**scenario)
                plots = generate_inline_plots(model)
                repeated = generate_inline_plots(model)
                self.assertEqual(plots, repeated)
                self.assertEqual(
                    {(str(item["run"]), str(item["kind"])) for item in plots},
                    expected_plot_set_from_fixture(model),
                )
                observed_run_counts.add(scenario["run_count"])
                observed_mark_counts.add(scenario["mark_count"])
                observed_barcode_counts.add(len(scenario["barcode_counts"]))
                retention_runs = {item["run"] for item in plots if item["kind"] == "retention"}
                self.assertEqual(retention_runs, {branch.run for branch in model.branches})
                if "DNA" in scenario["modalities"]:
                    complexity_runs = {item["run"] for item in plots if item["kind"] == "complexity"}
                    self.assertEqual(complexity_runs, {branch.run for branch in model.branches if branch.modality == "DNA"})
                else:
                    self.assertFalse(any(item["kind"] in {"complexity", "mark_composition"} for item in plots))
                expected_sample_categories = {
                    str(row["barcode_sequence"] or row["category"])
                    for row in model.barcode_composition
                    if row["barcode_type"] == "sample_barcode"
                }
                for plot in plots:
                    _, _, boxes, represented = svg_geometry(plot["svg"])
                    plot_width, plot_height, _, _ = svg_geometry(plot["svg"])
                    assert_non_overlapping_layout(self, boxes)
                    if plot["kind"] == "sample_composition":
                        sample_widths[len(scenario["barcode_counts"])] = plot_width
                        self.assertEqual(represented, expected_sample_categories)
                        rendered_text = "".join(ET.fromstring(plot["svg"]).itertext())
                        self.assertIn(str(plot["run"]), rendered_text)
                        root = ET.fromstring(plot["svg"])
                        statically_identified = {
                            element.attrib["data-static-category"]
                            for element in root.iter()
                            if "data-static-category" in element.attrib
                        }
                        categories_for_run = {
                            str(row["barcode_sequence"] or row["category"])
                            for row in model.barcode_composition
                            if row["run"] == plot["run"]
                            and row["barcode_type"] == "sample_barcode"
                        }
                        self.assertEqual(statically_identified, categories_for_run)
                    if plot["kind"] == "mark_composition":
                        expected_marks = {
                            str(row["category"])
                            for row in model.barcode_composition
                            if row["run"] == plot["run"]
                            and row["barcode_type"] == "dna_mark"
                            and row["category"] != "NoMatch"
                        }
                        self.assertEqual(represented, expected_marks)
                        self.assertNotIn("NoMatch", represented)
                        rendered_text = "".join(ET.fromstring(plot["svg"]).itertext())
                        self.assertIn(str(plot["run"]), rendered_text)
                        self.assertIn("Total:", rendered_text)
                        self.assertNotIn("denominator", rendered_text.lower())
                    if plot["kind"] == "retention":
                        retention_heights.setdefault(scenario["run_count"], plot_height)
        self.assertEqual(observed_run_counts, {1, 2, 4})
        self.assertEqual(observed_mark_counts, {1, 2, 6})
        self.assertEqual(observed_barcode_counts, {1, 4, 12, 30})
        self.assertEqual(sorted(sample_widths), [1, 4, 12, 30])
        self.assertEqual([sample_widths[count] for count in sorted(sample_widths)], sorted(sample_widths.values()))
        self.assertGreater(retention_heights[4], retention_heights[2])
        self.assertGreater(retention_heights[2], retention_heights[1])

    def test_sample_barcode_palette_scales_and_separates_adjacent_categories(self):
        labels = [f"barcode-{index:02d}" for index in range(30)] + ["NoMatch"]
        palette = categorical_palette(labels, "sample-barcode")
        self.assertEqual(palette, categorical_palette(reversed(labels), "sample-barcode"))
        self.assertEqual(palette["NoMatch"], "#aeb4b4")
        self.assertEqual(len(set(palette.values())), len(labels))
        hues = {
            label: float(re.match(r"hsl\(([0-9.]+)", color).group(1))
            for label, color in palette.items()
            if label != "NoMatch"
        }
        ordered = sorted(hues)
        for left, right in zip(ordered, ordered[1:]):
            distance = abs(hues[left] - hues[right]) % 360
            self.assertGreaterEqual(min(distance, 360 - distance), 7.9)

    def test_every_inline_svg_has_a_unique_download_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            model = synthetic_plot_model(
                run_count=4,
                modalities=("RNA", "DNA"),
                mark_count=6,
                barcode_counts=tuple([1] * 7 + [50] * 23),
                no_match_count=1_000,
                long_labels=True,
                roi=False,
            )
            plots = generate_inline_plots(model)
            report = Path(tmp) / "report.html"
            render_report(model, plots, report)
            document = report.read_text(encoding="utf-8")
            svg_ids = re.findall(r'<svg\s+id="([^"]+)"', document)
            controls = re.findall(r'<button[^>]+id="([^"]+)"[^>]+data-svg-download="([^"]+)"[^>]+data-filename="([^"]+)"', document)
            self.assertEqual(len(svg_ids), len(plots))
            self.assertEqual(len(svg_ids) * 2, len(controls))
            self.assertEqual(len(svg_ids), len(set(svg_ids)))
            self.assertEqual(len({control[0] for control in controls}), len(controls))
            self.assertEqual(len({control[2] for control in controls}), len(controls))
            self.assertEqual({control[1] for control in controls}, set(svg_ids))
            self.assertIn("new XMLSerializer()", document)
            self.assertIn("URL.createObjectURL(blob)", document)
            self.assertIn("canvas.toBlob", document)
            self.assertIn("clone.setAttribute('width',String(viewBox[2]))", document)
            self.assertIn("clone.setAttribute('height',String(viewBox[3]))", document)
            self.assertIn("context.fillRect(0,0,prepared.viewBox[2],prepared.viewBox[3])", document)
            self.assertIn("context.drawImage(image,0,0,prepared.viewBox[2],prepared.viewBox[3])", document)
            self.assertIn("Download SVG", document)
            self.assertIn("Download PNG", document)

    def test_overview_components_are_aggregated_from_generic_source_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            model = synthetic_plot_model(
                run_count=2,
                modalities=("RNA", "DNA"),
                mark_count=3,
                barcode_counts=(7, 5, 3, 1),
                no_match_count=2,
                roi=False,
            )
            report = Path(tmp) / "report.html"
            render_report(model, generate_inline_plots(model), report)
            document = report.read_text(encoding="utf-8")
            for run_index, run in enumerate(sorted({branch.run for branch in model.branches})):
                branch_rows = {}
                for row in model.library_complexity:
                    if row["run"] == run:
                        branch_rows.setdefault(row["branch"], row)
                examined = sum(int(row["read_pairs_examined"]) for row in branch_rows.values())
                pcr = sum(int(row["pcr_or_library_duplicate_pairs"]) for row in branch_rows.values())
                optical = sum(int(row["optical_duplicate_pairs"]) for row in branch_rows.values())
                self.assertIn(f"PCR/library: {pcr:,} ({pcr / examined * 100:.2f}% of examined pairs)", document)
                self.assertIn(f"Optical: {optical:,} ({optical / examined * 100:.2f}% of examined pairs)", document)
                saturation = 41.25 + run_index
                self.assertIn(f"{saturation:.1f}%", document)

    def test_missing_rna_saturation_omits_only_rna_complexity_plot(self):
        model = synthetic_plot_model(
            run_count=2,
            modalities=("RNA", "DNA"),
            mark_count=2,
            barcode_counts=(2, 2, 2, 2),
            no_match_count=1,
            roi=False,
            rna_saturation=False,
        )
        actual = {(str(item["run"]), str(item["kind"])) for item in generate_inline_plots(model)}
        self.assertEqual(actual, expected_plot_set_from_fixture(model))
        self.assertFalse(any(kind == "rna_complexity" for _, kind in actual))

    def test_rna_saturation_is_read_from_summary_without_unweighted_aggregation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=True, dna=False, groups={"run-a": ["BAR-A"]})
            model = build_report_model(root)
            saturation = [
                row for row in model.qc_metrics
                if row["metric"] == "sequencing_saturation_pct"
            ]
            self.assertEqual([row["value"] for row in saturation], [37.5])
            self.assertTrue(all("unique GeneFull feature" in str(row["denominator"]) for row in saturation))

            duplicate_source = dict(saturation[0])
            duplicate_source["branch"] = "RNA replicate"
            duplicate_source["source"] = "second.Summary.csv"
            duplicate_source["value"] = 52.0
            model.qc_metrics.append(duplicate_source)
            report = Path(tmp) / "report.html"
            render_report(model, generate_inline_plots(model), report)
            document = report.read_text(encoding="utf-8")
            self.assertIn("2 branch values; see RNA complexity", document)
            rna_plot = next(item for item in generate_inline_plots(model) if item["kind"] == "rna_complexity")
            rendered_text = "".join(ET.fromstring(rna_plot["svg"]).itertext())
            self.assertIn("37.50%", rendered_text)
            self.assertIn("52.00%", rendered_text)
            self.assertNotIn("44.75%", rendered_text)

    def test_invalid_rna_saturation_source_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            QcFixture(root, rna=True, dna=False, groups={"run-a": ["BAR-A"]})
            summary = next((root / "rna_align").glob("*.Solo.outGeneFull/Summary.csv"))
            write(summary, "Number of Reads,20\nSequencing Saturation,1.2\nEstimated Number of Cells,2\n")
            with self.assertRaisesRegex(ValueError, "saturation outside 0–1"):
                build_report_model(root)

    def test_exportable_svgs_have_explicit_background_and_format_filenames(self):
        with tempfile.TemporaryDirectory() as tmp:
            model = synthetic_plot_model(
                run_count=2,
                modalities=("RNA", "DNA"),
                mark_count=2,
                barcode_counts=(1, 3, 5, 7),
                no_match_count=11,
                roi=True,
            )
            renamed = {"run-1": "δοκιμή-με-μακρύ-όνομα", "run-2": "測試獨立執行名稱"}
            model.branches = [
                Branch(renamed[branch.run], branch.parent, branch.modality, branch.branch, branch.group, branch.split_id)
                for branch in model.branches
            ]
            for collection in (model.retention, model.qc_metrics, model.barcode_composition, model.library_complexity):
                for row in collection:
                    row["run"] = renamed[str(row["run"])]
            for row in model.barcode_composition:
                if row["barcode_type"] == "sample_barcode" and row["category"] != "NoMatch":
                    row["barcode_label"] = f"βιοδείκτης-με-μακρά-περιγραφή-{row['barcode_sequence']}"
            plots = generate_inline_plots(model)
            for item in plots:
                root = ET.fromstring(item["svg"])
                backgrounds = [
                    element for element in root
                    if element.attrib.get("data-export-background") == "white"
                ]
                self.assertEqual(len(backgrounds), 1)
            report = Path(tmp) / "report.html"
            render_report(model, plots, report)
            document = report.read_text(encoding="utf-8")
            filenames = re.findall(r'data-filename="([^"]+)"', document)
            self.assertEqual(len(filenames), len(plots) * 2)
            self.assertEqual(len(filenames), len(set(filenames)))
            stems = {}
            for name in filenames:
                stem, suffix = name.rsplit(".", 1)
                stems.setdefault(stem, set()).add(suffix)
            self.assertTrue(all(formats == {"svg", "png"} for formats in stems.values()))

    def test_roi_plot_has_clear_title_and_branch_legend(self):
        model = synthetic_plot_model(
            run_count=1,
            modalities=("DNA",),
            mark_count=3,
            barcode_counts=(5, 5),
            no_match_count=1,
            roi=True,
        )
        plot = next(
            item for item in generate_inline_plots(model)
            if item["kind"] == "complexity"
        )
        root = ET.fromstring(plot["svg"])
        rendered_text = "".join(root.itertext())

        self.assertIn(
            "Predicted unique pairs with deeper sequencing",
            rendered_text,
        )
        self.assertNotIn(
            "Picard ROI estimate (actual histogram points)",
            rendered_text,
        )

        legend_labels = sorted(
            element.attrib["data-roi-legend-label"]
            for element in root.iter()
            if "data-roi-legend-label" in element.attrib
        )
        expected = sorted(
            {
                branch.branch
                for branch in model.branches
                if branch.modality == "DNA"
            }
        )
        self.assertEqual(legend_labels, expected)

    def test_overview_uses_unspecified_when_dna_mode_metadata_is_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            model = synthetic_plot_model(
                run_count=1,
                modalities=("DNA",),
                mark_count=1,
                barcode_counts=(5,),
                no_match_count=1,
                roi=False,
            )
            self.assertEqual(model.run_metadata, {})

            report = Path(tmp) / "report.html"
            render_report(model, generate_inline_plots(model), report)
            document = report.read_text(encoding="utf-8")
            overview = document.split(
                '<section id="overview"', 1
            )[1].split("</section>", 1)[0]

            self.assertIn("parent-1: unspecified", overview)
            self.assertNotIn(
                "DNA tagmentation</span><strong>not present",
                overview,
            )

    def test_duplicate_explanation_matches_the_audited_pipeline_order(self):
        align_script = (REPO_ROOT / "scripts" / "core_runtime" / "AlignDNA.sh").read_text(encoding="utf-8")
        dna_workflow = (REPO_ROOT / "subworkflows" / "local" / "dna_core" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("--require-flags 0x2", align_script)
        self.assertIn("-L ${blacklist_bed}", align_script)
        self.assertLess(
            dna_workflow.index("GATK4_MARKDUPLICATES(ch_gatk_markduplicates_input"),
            dna_workflow.index("NORMALIZE_DNA_MARKDUPLICATES(ch_normalize_markduplicates_input)"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            model = synthetic_plot_model(
                run_count=1,
                modalities=("DNA",),
                mark_count=1,
                barcode_counts=(9,),
                no_match_count=1,
                roi=False,
            )
            report = Path(tmp) / "report.html"
            render_report(model, generate_inline_plots(model), report)
            document = report.read_text(encoding="utf-8")
            self.assertIn("Percent duplicates", document)
            self.assertIn("blacklist-overlap exclusion", document)
            self.assertIn("proper-pair flag requirement", document)
            self.assertIn("canonical-chromosome normalization occurs after MarkDuplicates", document)
            self.assertIn("no MAPQ threshold before this metric", document)

    def test_production_qc_code_contains_no_completed_run_identities(self):
        forbidden = (
            rf"\b{'group' + '1'}\b", rf"\b{'group' + '2'}\b",
            rf"\b{'AG' + 'CT'}\b", rf"\b{'G' + 'AC'}\b", rf"\b{'G' + 'AT'}\b",
            "H3" + "K27ac", "H3" + "K27me3", "H3" + "K9me3",
            "TreSeqDT_" + "260130_Tonsil",
        )
        production_paths = [
            *sorted((REPO_ROOT / "lib" / "tresflow_qc").glob("*.py")),
            REPO_ROOT / "bin" / "assess_tresflow_run.py",
            REPO_ROOT / "bin" / "render_tres_report.py",
            REPO_ROOT / "modules" / "local" / "tres_report_html" / "main.nf",
        ]
        for path in production_paths:
            source = path.read_text(encoding="utf-8")
            for pattern in forbidden:
                self.assertIsNone(re.search(pattern, source), f"{pattern} found in {path}")


if __name__ == "__main__":
    unittest.main()
