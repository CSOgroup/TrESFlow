#!/usr/bin/env python3

"""Resolve human canonical chromosomes without renaming reference contigs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


AUTOSOMES = tuple(str(number) for number in range(1, 23))
UCSC_NAMES = frozenset(
    [*(f"chr{chromosome}" for chromosome in AUTOSOMES), "chrX", "chrY", "chrM"]
)
ENSEMBL_NAMES = frozenset([*AUTOSOMES, "X", "Y", "MT", "M"])


class CanonicalChromosomeError(ValueError):
    """Raised when a reference dictionary is not safely classifiable."""


def read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    entries: list[tuple[str, int]] = []
    seen: set[str] = set()
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise CanonicalChromosomeError(
                f"Malformed chromosome-size line {line_number} in {path}: {raw_line!r}"
            )
        name = fields[0]
        try:
            length = int(fields[1])
        except ValueError as error:
            raise CanonicalChromosomeError(
                f"Invalid chromosome length on line {line_number} in {path}: {fields[1]!r}"
            ) from error
        if length < 1:
            raise CanonicalChromosomeError(
                f"Chromosome length must be positive on line {line_number} in {path}: {length}"
            )
        if name in seen:
            raise CanonicalChromosomeError(f"Duplicate chromosome {name!r} in {path}")
        seen.add(name)
        entries.append((name, length))
    if not entries:
        raise CanonicalChromosomeError(f"Reference chromosome dictionary is empty: {path}")
    return entries


def read_bwa_ann(path: Path) -> list[tuple[str, int]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise CanonicalChromosomeError(f"BWA annotation file is empty: {path}")
    header = lines[0].split()
    try:
        sequence_count = int(header[1])
    except (IndexError, ValueError) as error:
        raise CanonicalChromosomeError(
            f"Cannot read the sequence count from BWA annotation header in {path}: {lines[0]!r}"
        ) from error
    expected_lines = 1 + (sequence_count * 2)
    if sequence_count < 1 or len(lines) < expected_lines:
        raise CanonicalChromosomeError(
            f"Malformed BWA annotation file {path}: expected {sequence_count} sequence entries"
        )

    entries: list[tuple[str, int]] = []
    seen: set[str] = set()
    for index in range(sequence_count):
        descriptor = lines[1 + (index * 2)].split()
        coordinates = lines[2 + (index * 2)].split()
        try:
            name = descriptor[1]
            length = int(coordinates[1])
        except (IndexError, ValueError) as error:
            raise CanonicalChromosomeError(
                f"Malformed BWA sequence entry {index + 1} in {path}"
            ) from error
        if length < 1:
            raise CanonicalChromosomeError(
                f"BWA chromosome length must be positive for {name!r} in {path}: {length}"
            )
        if name in seen:
            raise CanonicalChromosomeError(f"Duplicate chromosome {name!r} in {path}")
        seen.add(name)
        entries.append((name, length))
    return entries


def resolve_canonical_entries(
    entries: list[tuple[str, int]], source_label: str
) -> tuple[str, list[tuple[str, int]]]:
    names = {name for name, _ in entries}
    ucsc_hits = names & UCSC_NAMES
    ensembl_hits = names & ENSEMBL_NAMES

    if ucsc_hits and ensembl_hits:
        raise CanonicalChromosomeError(
            "Cannot determine a coherent human chromosome naming convention for "
            f"{source_label}: both UCSC-style ({', '.join(sorted(ucsc_hits))}) and "
            f"Ensembl-style ({', '.join(sorted(ensembl_hits))}) canonical names are present"
        )
    if not ucsc_hits and not ensembl_hits:
        raise CanonicalChromosomeError(
            "Cannot determine a human chromosome naming convention for "
            f"{source_label}: no exact canonical autosome, X, Y, or mitochondrial names were found"
        )

    if ucsc_hits:
        style = "ucsc"
        allowed_names = UCSC_NAMES
        missing_anchors = [
            label
            for label, present in (
                ("an autosome (chr1-chr22)", any(f"chr{chrom}" in names for chrom in AUTOSOMES)),
                ("chrX", "chrX" in names),
                ("chrY", "chrY" in names),
                ("chrM", "chrM" in names),
            )
            if not present
        ]
    else:
        style = "ensembl"
        mitochondrial = [name for name in ("MT", "M") if name in names]
        if len(mitochondrial) > 1:
            raise CanonicalChromosomeError(
                f"Cannot choose a single Ensembl mitochondrial chromosome for {source_label}: "
                "both MT and M are present"
            )
        allowed_names = frozenset([*AUTOSOMES, "X", "Y", *mitochondrial])
        missing_anchors = [
            label
            for label, present in (
                ("an autosome (1-22)", any(chrom in names for chrom in AUTOSOMES)),
                ("X", "X" in names),
                ("Y", "Y" in names),
                ("MT or M", bool(mitochondrial)),
            )
            if not present
        ]

    if missing_anchors:
        raise CanonicalChromosomeError(
            f"Cannot safely resolve canonical human chromosomes for {source_label}; missing: "
            + ", ".join(missing_anchors)
        )

    canonical = [(name, length) for name, length in entries if name in allowed_names]
    return style, canonical


def verify_matching_contracts(
    primary: tuple[str, list[tuple[str, int]]],
    secondary: tuple[str, list[tuple[str, int]]],
    primary_label: str,
    secondary_label: str,
) -> None:
    primary_style, primary_entries = primary
    secondary_style, secondary_entries = secondary
    if primary_style != secondary_style:
        raise CanonicalChromosomeError(
            f"Reference naming mismatch: {primary_label} is {primary_style}-style but "
            f"{secondary_label} is {secondary_style}-style"
        )
    if dict(primary_entries) != dict(secondary_entries):
        raise CanonicalChromosomeError(
            f"Canonical chromosome names or lengths disagree between {primary_label} and "
            f"{secondary_label}"
        )


def write_contract(
    modality: str,
    source: Path,
    source_format: str,
    output_dir: Path,
    verification_sizes: Path | None = None,
) -> dict[str, object]:
    entries = read_bwa_ann(source) if source_format == "bwa-ann" else read_chrom_sizes(source)
    resolved = resolve_canonical_entries(entries, str(source))
    if verification_sizes:
        verification_entries = read_chrom_sizes(verification_sizes)
        verification = resolve_canonical_entries(verification_entries, str(verification_sizes))
        verify_matching_contracts(resolved, verification, str(source), str(verification_sizes))

    style, canonical_entries = resolved
    output_dir.mkdir(parents=True, exist_ok=True)
    allowlist = output_dir / f"{modality}_canonical_chromosomes.txt"
    chrom_sizes = output_dir / f"{modality}_canonical_chromosomes.chrom.sizes"
    allowlist.write_text(
        "".join(f"{name}\n" for name, _ in canonical_entries), encoding="utf-8"
    )
    chrom_sizes.write_text(
        "".join(f"{name}\t{length}\n" for name, length in canonical_entries),
        encoding="utf-8",
    )
    return {
        "style": style,
        "source": str(source.resolve()),
        "allowlist": str(allowlist.resolve()),
        "chrom_sizes": str(chrom_sizes.resolve()),
        "contigs": [name for name, _ in canonical_entries],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--rna-chrom-sizes", type=Path)
    parser.add_argument("--dna-bwa-ann", type=Path)
    parser.add_argument("--dna-chrom-sizes", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.rna_chrom_sizes and not args.dna_bwa_ann:
        raise CanonicalChromosomeError("At least one RNA or DNA reference dictionary is required")

    contracts: dict[str, object] = {}
    if args.rna_chrom_sizes:
        contracts["rna"] = write_contract(
            "rna", args.rna_chrom_sizes, "chrom-sizes", args.output_dir
        )
    if args.dna_bwa_ann:
        contracts["dna"] = write_contract(
            "dna",
            args.dna_bwa_ann,
            "bwa-ann",
            args.output_dir,
            args.dna_chrom_sizes,
        )
    print(json.dumps(contracts, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CanonicalChromosomeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
