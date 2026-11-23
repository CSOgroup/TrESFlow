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
    assert "seq11" in html.split("Exact barcode counts and denominators")[1]
    assert html.count("data-svg-download=") == 2


def test_concise_cards_keep_counts_and_correct_denominators():
    from bulk_report import target_card, saturation
    rna = dict(routed_pairs=50, molecules=10, featureU=40, genes_detected=7)
    assert saturation(rna) == 75
    html = target_card("RNA", rna, True, raw=100)
    assert "10 (10.00% of raw modality pairs)" in html
    assert "75.00%" in html and "Pools with molecules" not in html
    supplied = target_card("RNA", rna, True)
    assert "10 (20.00% of supplied group pairs)" in supplied
    assert "unavailable" in supplied and "Supplied group pairs" in supplied
    assert saturation(dict(featureU=0, molecules=0)) is None
    dna = dict(usable_fragments=80, bwa_primary_pairs=100, read_pair_duplicates=25,
               read_pairs_examined=100, peak_count=2, fragments_in_peaks=20)
    html = target_card("DNA", dna, False, raw=200)
    assert "80 (40.00% of raw modality pairs)" in html
    assert "FRiP · final aligned fragments" in html and "25.00%" in html


def test_picard_prediction_uses_actual_histogram(tmp_path):
    from bulk_report import dna_complexity
    q = dict(picard=dict(read_pairs_examined=100, read_pair_duplicates=30,
                        read_pair_optical_duplicates=5, estimated_library_size=200,
                        roi=[dict(depth_multiplier=1, yield_multiplier=1),
                             dict(depth_multiplier=2, yield_multiplier=1.5)]))
    html = dna_complexity(q, tmp_path, "id", "sample", "mark")
    assert "Observed duplicate composition" in html
    assert "mark 2×: 105 estimated unique pairs" in html
    assert "PCR/library 25 (25.00%)" in html
    q["picard"]["roi"] = []
    assert "Picard supplied no ROI histogram" in dna_complexity(q, tmp_path, "id", "sample", "mark")
    q["picard"]["read_pairs_examined"] = 0
    assert "composition and depth prediction unavailable" in dna_complexity(q, tmp_path, "id", "sample", "mark")


def test_ligation_plots_are_collapsed():
    panel = BarcodePanel("s · RNA", "rna", "L1", "Split input", 3, "split-input pairs",
                         Counter({("AAA", ""): 3}))
    html = barcode_section([panel])
    assert html.startswith('<details><summary>Ligation barcode abundance')
    assert '<details open' not in html


def test_fragment_distribution_preserves_denominator_and_tail():
    from tresflow_qc.bulk_fragments import fragment_figure
    svg = fragment_figure({"100": 6, "200": 2, "2000": 2}, "sample · DNA")["svg"]
    ElementTree.fromstring(svg)
    assert "6 fragments (60.00%)" in svg
    assert "20.00% beyond 1,000 bp" in svg
    assert "Median 100 bp" in svg and "5 bp bins" in svg
    full = fragment_figure({"100": 6, "200": 2, "2000": 2}, "sample · DNA", True)["svg"]
    assert "2 fragments (20.00%)" in full and "0.00% beyond 2,000 bp" in full
    ElementTree.fromstring(fragment_figure({}, "empty")["svg"])
    assert "Median 150 bp" in fragment_figure({"100": 1, "200": 1}, "even")["svg"]


def test_raw_sample_plot_and_percentage_warnings(tmp_path):
    from bulk_report import preprocessing
    prefix = tmp_path / "sample.rna"
    Path(str(prefix) + "_sample_barcode.stats.tsv").write_text(
        "reads\t100\tinput pairs\nreads_without_bc\t10\tinput pairs\n")
    Path(str(prefix) + "_sample_barcode.counts.tsv").write_text("90\tATGC\n10\tNoMatch\n")
    Path(str(prefix) + "_barcode_gates.tsv").write_text("metric\tpairs\nligation_accepted_pairs\t50\n")
    Path(str(prefix) + "_barcode_composition.tsv").write_text(
        "sample_id\tmodality\tbarcode_type\tgroup\tcategory\tcount\tdenominator_count\tdenominator_definition\n"
        "sample\trna\tsample_barcode\tgroup\tATGC\t40\t50\tligation-accepted pairs\n"
        "sample\trna\tsample_barcode\t__all__\tNoMatch\t10\t50\tligation-accepted pairs\n")
    panels, tables, warnings = preprocessing(tmp_path)
    assert len(panels) == 1 and panels[0].stage == "Sample-barcode tagging"
    assert panels[0].denominator == 100
    assert "10.00% of raw input pairs" in warnings[0]
    assert "20.00% of ligation-accepted pairs" in warnings[1]
    assert "10 of" not in " ".join(warnings)
    assert "Raw input pairs" in tables[0]


def test_final_fragments_include_pairs_outside_peaks(tmp_path, monkeypatch):
    import pysam
    from bulk_dna import create_fragments, overlap_count
    monkeypatch.chdir(tmp_path)
    header = dict(HD=dict(VN="1.6"), SQ=[dict(SN="chr1", LN=1000)], RG=[dict(ID="rg")])
    with pysam.AlignmentFile("pairs.bam", "wb", header=header) as bam:
        for name, start, duplicate in (("a", 100, 0), ("b", 300, 0x400)):
            for flag, pos, mate in ((99, start, start+20), (147, start+20, start)):
                read = pysam.AlignedSegment(bam.header)
                read.query_name = name
                read.flag = flag | duplicate
                read.reference_id = read.next_reference_id = 0
                read.reference_start, read.next_reference_start = pos, mate
                read.cigarstring = "20M"
                read.mapping_quality = 60
                for tag, value in (("RG", "rg"), ("CB", "cell"), ("SB", "sample")):
                    read.set_tag(tag, value)
                bam.write(read)
    lengths, duplicates = create_fragments("pairs.bam", {"chr1": 1000})
    Path("peaks.bed").write_text("chr1\t100\t140\n")
    assert sum(lengths.values()) == 2 and duplicates == 1
    assert overlap_count(Path("peaks.bed")) == (1, 40)
    with pysam.AlignmentFile("spans.bam") as bam:
        assert sum(1 for _ in bam) == 2
