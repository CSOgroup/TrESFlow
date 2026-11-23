"""Read-pair retention stages observable in raw and routed bulk analysis."""
from itertools import groupby
from pathlib import Path

from .core import _metric_lookup, _retention_lookup, read_stats, read_tsv
from .plots import DNA, RNA, mark_palette, retention_figure


def stage_rows(stages, denominator, basis, scope):
    rows = []
    previous = None
    for stage, count, source in stages:
        if count is None:
            continue
        count = int(count)
        if count < 0 or (previous is not None and count > previous):
            raise ValueError(f"{scope}: non-nested retention at {stage} ({count}), source {source}")
        rows.append(dict(stage=stage, stage_order=len(rows), output_pairs=count, count_source=source,
                         cumulative_raw_pct=100 * count / denominator if denominator else None,
                         percentage_basis=basis if denominator else None))
        previous = count
    return rows


def shared_stages(root, sample, modality, mode):
    prefix = f"{sample}.{modality}"
    stats_path = root / f"{prefix}_sample_barcode.stats.tsv"
    raw = read_stats(stats_path).get("reads", {}).get("count") if stats_path.exists() else None
    stages = [("Raw input", raw, f"{stats_path.name}:reads")]
    split_path = root / f"{prefix}_read_retention.tsv"
    split = _retention_lookup(read_tsv(split_path)) if split_path.exists() else {}
    gate_path = root / f"{prefix}_barcode_gates.tsv"
    gates = _metric_lookup(read_tsv(gate_path)) if gate_path.exists() else {}
    split_input = split.get(("__all__", "__all__", "split_input_pairs"))
    if split_input is not None and gates.get("split_input_pairs", split_input) != split_input:
        raise ValueError(f"{split_path}: split input disagrees with {gate_path.name}")
    entering = gates.get("split_input_pairs", split_input)
    artifact_path = root / f"{sample}.dual_tag_artifact_filter.summary.tsv"
    if modality == "dna" and artifact_path.exists():
        artifact, = read_tsv(artifact_path)
        kept = int(artifact["retained_pairs"])
        if entering is not None and kept != entering:
            raise ValueError(f"{artifact_path}: retained pairs disagree with split input")
        stages.extend([
            ("After paired trimming", int(artifact["input_pairs"]), f"{artifact_path.name}:input_pairs"),
            ("After DT artifact filter", kept, f"{artifact_path.name}:retained_pairs"),
        ])
    else:
        label = "After paired trimming" if modality == "rna" or mode == "single" else "Pairs entering routing"
        stages.append((label, entering, f"{gate_path.name if gates else split_path.name}:split_input_pairs"))
    for label, field in (
        ("L1–L2–L3 accepted", "ligation_barcode_accepted_pairs"),
        ("Sample barcode accepted", "sample_barcode_accepted_pairs"),
        ("Mark barcode accepted", "modality_barcode_accepted_pairs"),
    ):
        if field in gates:
            stages.append((label, gates[field], f"{gate_path.name}:{field}"))
    stages = [row for row in stages if row[1] is not None]
    return stages, raw, split


def retention_figures(targets, observations, root, restart):
    """Show shared stages once per sample/modality and each downstream group."""
    root = Path(root)
    colors = mark_palette(t["mark"] for t in targets if t["modality"] == "dna")
    figures = []
    for sample, sample_targets in groupby(sorted(targets, key=lambda t: t["sample"]), key=lambda t: t["sample"]):
        panels = []
        for modality, subset in groupby(sorted(sample_targets, key=lambda t: t["modality"] != "rna"),
                                       key=lambda t: t["modality"]):
            subset = sorted(subset, key=lambda t: (t["group"], t.get("mark") or ""))
            shared, raw, split = ([], None, {}) if restart else shared_stages(
                root, sample, modality, subset[0].get("dna_tagmentation"))
            basis = f"raw {modality.upper()} input"
            if shared:
                panels.append((stage_rows(shared, raw, basis, sample), RNA if modality == "rna" else DNA,
                               f"{modality.upper()} shared processing · % of {basis}" if raw else
                               f"{modality.upper()} shared processing · percentages unavailable"))
            for target in subset:
                q = observations.get((modality, target["id"]))
                if q is None:
                    continue
                source = f"{target['id']}.{modality}.qc.json"
                rna = modality == "rna"
                first = q["routed_pairs"] if rna else q["bwa_primary_pairs"]
                denominator = raw if raw is not None else first
                measured = "supplied" if restart else "routed"
                basis = f"raw {modality.upper()} input" if raw is not None else (
                    f"{measured} group pairs" if rna else f"{measured} mark pairs")
                stages = shared[-1:] if shared else []
                group_count = None
                if not restart:
                    group_count = split.get((target["group"], "__all__", "routed_group_pairs"))
                    stages = stages + [("Group routing", group_count, f"{sample}.{modality}_read_retention.tsv:routed_group_pairs")]
                if rna:
                    if group_count != first:
                        stages.append(("Routed group pairs", first, f"{source}:routed_pairs"))
                    stages += [
                        ("STAR input", q["retained_pairs"], f"{source}:retained_pairs"),
                        ("Mapped pairs", q["genomeU"] + q["genomeM"], f"{source}:genomeU+genomeM"),
                        ("Gene-assigned pairs", q["featureU"], f"{source}:featureU"),
                    ]
                else:
                    branch_count = split.get((target["group"], target["mark"], "routed_branch_pairs"))
                    if not restart and branch_count is not None and branch_count != first:
                        stages.append(("Mark routing", branch_count, f"{sample}.dna_read_retention.tsv:routed_branch_pairs"))
                    first_label = "Mark routing" if not restart and branch_count == first else "Alignment input"
                    stages += [
                        (first_label, first, f"{source}:bwa_primary_pairs"),
                        ("Blacklist filter", q["post_blacklist_primary_pairs"], f"{source}:post_blacklist_primary_pairs"),
                        ("Mapped pairs (R1)", q.get("post_blacklist_mapped_primary_pairs"), f"{source}:post_blacklist_mapped_primary_pairs"),
                        ("Proper pairs", q["proper_pair_primary_pairs"], f"{source}:proper_pair_primary_pairs"),
                        ("Final fragments", q["usable_fragments"], f"{source}:usable_fragments"),
                    ]
                title = " · ".join([modality.upper(), target["group"]] + ([] if rna else [target["mark"]]))
                title += f" · % of {basis}" if denominator else " · percentages unavailable"
                panels.append((stage_rows(stages, denominator, basis, target["id"]),
                               RNA if rna else colors[target["mark"]], title))
        if panels:
            figures.append(retention_figure(sample, panels))
    return figures
