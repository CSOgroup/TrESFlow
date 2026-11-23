"""Bulk retention uses measurable pair counts with the shared report renderer."""
import sys
from pathlib import Path
from xml.etree import ElementTree

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from bulk_report import run_overview
from tresflow_qc.bulk_retention import retention_figures, shared_stages, stage_rows


def test_overview_has_run_scope_and_available_metadata():
    targets = [dict(sample="sample", group="g1", modality="rna"),
               dict(sample="sample", group="g1", modality="dna", mark="H3K27ac", dna_tagmentation="single"),
               dict(sample="sample", group="g1", modality="dna", mark="H3K27me3", dna_tagmentation="single")]
    html = run_overview(dict(pipeline_version="v1.2.0-test", generated_at="2026-11-23T12:00:00+00:00"), targets)
    root = ElementTree.fromstring(html)
    fields = {cell.find("span").text: cell.find("strong").text for cell in root}
    assert fields == {"Independent runs": "g1", "Modalities": "DNA, RNA",
                      "DNA marks": "H3K27ac, H3K27me3", "DNA tagmentation": "sample: single",
                      "Pipeline version": "v1.2.0-test", "Generated": "2026-11-23T12:00:00+00:00"}
    targets[-1] = dict(targets[-1], sample="other", dna_tagmentation=None)
    unknown = run_overview({}, targets)
    assert "sample_g1" in unknown and "other_g1" in unknown
    assert "other: unavailable" in unknown and "<strong>unavailable</strong>" in unknown


def test_pair_stages_keep_units_and_validate_subset_counts():
    rows = stage_rows([("Input", 100, "input"), ("Missing", None, ""), ("Assigned", 30, "gene")],
                      200, "raw RNA input", "sample")
    assert [r["output_pairs"] for r in rows] == [100, 30]
    assert rows[-1]["cumulative_raw_pct"] == 15
    assert rows[-1]["percentage_basis"] == "raw RNA input"
    empty = stage_rows([("Input", 0, "input")], 0, "raw", "sample")
    assert empty[0]["percentage_basis"] is None and empty[0]["cumulative_raw_pct"] is None
    from tresflow_qc.plots import retention_figure
    figure = retention_figure("zero", [(empty, "#000000", "RNA")])
    assert "Input: 0 pairs; percentage unavailable" in figure["svg"]
    assert "of None" not in figure["svg"]
    with pytest.raises(ValueError, match="non-nested retention"):
        stage_rows([("Input", 10, "input"), ("Assigned", 11, "gene")], 10, "raw", "sample")


def write_shared(root, modality, *, dual=False):
    prefix = f"s.{modality}"
    (root / f"{prefix}_sample_barcode.stats.tsv").write_text("reads\t100\t100\n")
    values = dict(split_input_pairs=80 if dual else 90, ligation_barcode_accepted_pairs=70,
                  sample_barcode_accepted_pairs=60)
    if modality == "dna":
        values["modality_barcode_accepted_pairs"] = 50
    (root / f"{prefix}_barcode_gates.tsv").write_text("metric\tpairs\n" +
        "".join(f"{key}\t{value}\n" for key, value in values.items()))
    split = "group\tbranch\tmetric\tpairs\n__all__\t__all__\tsplit_input_pairs\t" + str(values["split_input_pairs"]) + "\n"
    split += "g1\t__all__\trouted_group_pairs\t" + ("30" if modality == "dna" else "40") + "\n"
    split += "g2\t__all__\trouted_group_pairs\t20\n"
    if modality == "dna":
        split += "g1\tA\trouted_branch_pairs\t20\ng1\tB\trouted_branch_pairs\t10\ng2\tA\trouted_branch_pairs\t20\n"
    (root / f"{prefix}_read_retention.tsv").write_text(split)
    if dual:
        (root / "s.dual_tag_artifact_filter.summary.tsv").write_text(
            "input_pairs\tretained_pairs\n90\t80\n")


def test_raw_and_fastq_retention_keep_all_groups_and_duplicate_fragments(tmp_path):
    write_shared(tmp_path, "rna")
    write_shared(tmp_path, "dna", dual=True)
    targets, qcs = [], {}
    for group, count in (("g1", 40), ("g2", 20)):
        name = f"s_{group}"
        targets.append(dict(sample="s", group=group, modality="rna", id=name))
        qcs["rna", name] = dict(routed_pairs=count, retained_pairs=count-2, genomeU=count-7,
                               genomeM=2, featureU=count-8, molecules=3)
    for group, mark, count in (("g1", "A", 20), ("g1", "B", 10), ("g2", "A", 20)):
        name = f"s_{group}_{mark}"
        targets.append(dict(sample="s", group=group, modality="dna", id=name, mark=mark,
                            dna_tagmentation="dual"))
        qcs["dna", name] = dict(bwa_primary_pairs=count, post_blacklist_primary_pairs=count-1,
                               post_blacklist_mapped_primary_pairs=count-2, proper_pair_primary_pairs=count-3,
                               usable_fragments=count-4, usable_duplicate_fragments=count-5)
    figure, = retention_figures(targets, qcs, tmp_path, False)
    ElementTree.fromstring(figure["svg"])
    assert figure["svg"].count('class="panel-title"') == 7
    assert "RNA · g1" in figure["svg"] and "RNA · g2" in figure["svg"]
    assert "DNA · g1 · A" in figure["svg"] and "DNA · g2 · A" in figure["svg"]
    assert "After DT artifact filter: 80 pairs" in figure["svg"]
    assert "Gene-assigned pairs: 32 pairs; 32.0% of raw RNA input" in figure["svg"]
    assert "Final fragments: 16 pairs; 16.0% of raw DNA input" in figure["svg"]
    assert "molecules" not in figure["svg"]
    supplied, = retention_figures(targets, qcs, tmp_path, True)
    assert "Raw input" not in supplied["svg"] and "shared processing" not in supplied["svg"]
    assert "Gene-assigned pairs: 32 pairs; 80.0% of supplied group pairs" in supplied["svg"]
    assert "Final fragments: 16 pairs; 80.0% of supplied mark pairs" in supplied["svg"]
    assert "of raw" not in supplied["svg"]


def test_dual_stages_do_not_invent_trimming_counts_and_check_sources(tmp_path):
    write_shared(tmp_path, "dna", dual=True)
    summary = tmp_path / "s.dual_tag_artifact_filter.summary.tsv"
    summary.unlink()
    stages, raw, _ = shared_stages(tmp_path, "s", "dna", "dual")
    assert raw == 100 and stages[1][:2] == ("Pairs entering routing", 80)
    summary.write_text("input_pairs\tretained_pairs\n90\t79\n")
    with pytest.raises(ValueError, match="retained pairs disagree"):
        shared_stages(tmp_path, "s", "dna", "dual")
