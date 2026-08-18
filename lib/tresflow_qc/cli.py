"""Command-line entry point shared by pipeline and standalone assessment."""

from __future__ import annotations

import argparse
import base64
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

from .core import build_report_model, write_model_json, write_model_tables
from .html_report import render_report
from .plots import generate_inline_plots, write_inline_plots


FIGURE_SUFFIXES = (
    "_read_retention.png", "_read_retention.svg",
    "_sample_barcode_composition.png", "_sample_barcode_composition.svg",
    "_dna_mark_composition.png", "_dna_mark_composition.svg",
    "_dna_library_complexity.png", "_dna_library_complexity.svg",
    "_rna_library_complexity.png", "_rna_library_complexity.svg",
    "sample_barcode_composition_rna.svg", "sample_barcode_composition_dna.svg",
    "dna_mark_composition.svg",
)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build exact TrESFlow QC tables and a self-contained offline HTML report.")
    parser.add_argument("tresflow_output", nargs="?", type=Path, help="completed output tree or staged explicit report inputs")
    parser.add_argument("--input-dir", type=Path, help="alias for the positional input path (used by Nextflow)")
    parser.add_argument("--output-dir", type=Path, help="destination directory (default: outputs/<input-name>)")
    parser.add_argument("--output-html", default="tres_report.html", help="HTML filename or path")
    parser.add_argument("--output-json", default="tres_report_metrics.json", help="normalized JSON filename or path")
    parser.add_argument(
        "--title", "--library-name", dest="title", default="",
        help="report title override (default: basename of the assessed pipeline-output directory)",
    )
    parser.add_argument("--pipeline-version", default="unknown", help="TrESFlow version recorded in provenance")
    parser.add_argument("--run-metadata-json", help="JSON object with additional run metadata")
    parser.add_argument("--run-metadata-base64", help="base64-encoded JSON object (safe for workflow command transport)")
    parser.add_argument("--no-json", action="store_true", help="do not write the optional standalone normalized JSON")
    parser.add_argument("--no-standalone-figures", action="store_true", help="embed SVG figures in HTML without writing separate SVG files")
    parser.add_argument("--no-plots", action="store_true", help="write normalized tables/optional JSON without figures or HTML")
    return parser.parse_args(argv)


def _resolve_output(output_dir: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else output_dir / path


def _print_console(model, paths: list[Path]) -> None:
    grouped = defaultdict(list)
    for row in model.retention:
        grouped[(row["run"], row["modality"], row["branch"])].append(row)
    for (run, modality, branch), rows in sorted(grouped.items()):
        final = max(rows, key=lambda row: int(row["stage_order"]))
        print(f"{run} {modality} {branch}: {int(final['output_pairs']):,} pairs ({float(final['cumulative_raw_pct']):.1f}% of raw)")
    for warning in model.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    print("Outputs:")
    for path in paths:
        print(f"  {path.resolve()}")


def _clear_generated_figures(output_dir: Path) -> None:
    """Remove only stale assessor-owned figures before a deterministic rerun."""
    for path in output_dir.iterdir():
        if path.is_file() and path.name.endswith(FIGURE_SUFFIXES):
            path.unlink()


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    if args.tresflow_output and args.input_dir:
        print("error: provide either positional tresflow_output or --input-dir, not both", file=sys.stderr)
        return 2
    input_root = (args.input_dir or args.tresflow_output)
    if input_root is None:
        print("error: an input directory is required", file=sys.stderr)
        return 2
    input_root = input_root.expanduser().resolve()
    if not input_root.is_dir():
        print(f"error: not a directory: {input_root}", file=sys.stderr)
        return 2
    output_dir = (args.output_dir or (Path.cwd() / "outputs" / input_root.name)).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {}
    metadata_text = args.run_metadata_json
    if args.run_metadata_base64:
        if metadata_text:
            print("error: provide only one run-metadata encoding", file=sys.stderr)
            return 2
        try:
            metadata_text = base64.b64decode(args.run_metadata_base64).decode("utf-8")
        except Exception as error:
            print(f"error: invalid --run-metadata-base64: {error}", file=sys.stderr)
            return 2
    if metadata_text:
        try:
            metadata = json.loads(metadata_text)
        except json.JSONDecodeError as error:
            print(f"error: invalid --run-metadata-json: {error}", file=sys.stderr)
            return 2
        if not isinstance(metadata, dict):
            print("error: --run-metadata-json must be a JSON object", file=sys.stderr)
            return 2
    model = build_report_model(
        input_root,
        library_name=args.title or input_root.name,
        pipeline_version=args.pipeline_version,
        run_metadata=metadata,
    )
    if not model.retention:
        print("error: no valid sequential retention paths could be reconstructed", file=sys.stderr)
        return 1
    tables = write_model_tables(model, output_dir)
    json_path = None
    if not args.no_json:
        json_path = _resolve_output(output_dir, args.output_json)
        write_model_json(model, json_path)
    if not args.no_plots:
        _clear_generated_figures(output_dir)
    plots = [] if args.no_plots else generate_inline_plots(model)
    standalone_plots = []
    if plots and not args.no_standalone_figures:
        standalone_plots = write_inline_plots(plots, output_dir)
    html_path = _resolve_output(output_dir, args.output_html)
    if not args.no_plots:
        render_report(model, plots, html_path)
    outputs = list(tables.values())
    if json_path is not None:
        outputs.append(json_path)
    outputs.extend(Path(item["svg_path"]) for item in standalone_plots)
    if not args.no_plots:
        outputs.append(html_path)
    _print_console(model, outputs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
