#!/usr/bin/env python3

import argparse
import os
import shutil
from pathlib import Path


def ensure_clean_path(path: Path):
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def copy_directory(src: Path, dst: Path):
    ensure_clean_path(dst)
    shutil.copytree(src, dst, symlinks=False)


def link_or_copy_file(src: Path, dst: Path):
    ensure_clean_path(dst)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def infer_rna_samples(rna_solo_dirs, rna_filtered_bams):
    out = set()

    for solo_dir in rna_solo_dirs:
        name = solo_dir.name
        suffix = ".Solo.outGeneFull"
        if name.endswith(suffix):
            out.add(name[: -len(suffix)])

    for bam in rna_filtered_bams:
        name = bam.name
        suffix = ".filtered_cells.bam"
        if name.endswith(suffix):
            out.add(name[: -len(suffix)])

    return sorted(out)


def infer_dna_group_keys(dna_nodup_bams):
    out = set()

    for bam in dna_nodup_bams:
        name = bam.name
        suffix = "_NoDup.bam"
        if not name.endswith(suffix):
            continue
        stem = name[: -len(suffix)]
        tokens = stem.split("_")
        if len(tokens) < 3:
            continue
        out.add("_".join(tokens[:-1]))

    return sorted(out)


def write_pairs(path: Path, dna_group_keys, rna_samples):
    rna_sample_set = set(rna_samples)
    exact_pairs = [(key, key) for key in dna_group_keys if key in rna_sample_set]

    with open(path, "wt", encoding="utf-8") as handle:
        for dna_group_key, rna_sample in exact_pairs:
            handle.write(f"{dna_group_key}\t{rna_sample}\n")

    return exact_pairs


def write_manifest(path: Path, rows):
    with open(path, "wt", encoding="utf-8") as handle:
        handle.write("kind\tsource_path\tstaged_name\n")
        for kind, source, staged_name in rows:
            handle.write(f"{kind}\t{source}\t{staged_name}\n")


def write_readiness_note(path: Path, species: str, genome: str, exact_pairs):
    lines = [
        "Shared sc_process staging boundary",
        "",
        f"species\t{species}",
        f"genome\t{genome}",
        f"exact_pairs\t{len(exact_pairs)}",
        "",
        "Current readiness:",
        "- DNA NoDup BAMs are staged with launcher-style names.",
        "- RNA STARsolo directories are staged with launcher-style names.",
        "- pairs.tsv contains exact dna_group_key -> rna_sample mappings when grouped names match.",
        "",
        "Current blocker before one shared sc_process.py call on this server:",
        "- sc_process.py triggers SnapATAC2 TSSE annotation loading.",
        "- Without a local cached/provided gene annotation, SnapATAC2 attempts a network download.",
        "- The last observed missing remote file for the human/hg38 path was:",
        "  gencode.v41.basic.annotation.gff3.gz from ftp.ebi.ac.uk",
    ]

    with open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Stage flat shared inputs for one future sc_process.py call."
    )
    parser.add_argument("--stage-dir", required=True, type=Path)
    parser.add_argument("--mo-map", required=True, type=Path)
    parser.add_argument("--sb-group-map", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--rna-solo-dir", action="append", default=[], type=Path)
    parser.add_argument("--rna-filtered-bam", action="append", default=[], type=Path)
    parser.add_argument("--dna-nodup-bam", action="append", default=[], type=Path)
    args = parser.parse_args()

    stage_dir = args.stage_dir.resolve()
    stage_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []

    mo_map_dst = stage_dir / "mo_map.tsv"
    shutil.copy2(args.mo_map, mo_map_dst)
    manifest_rows.append(("map", str(args.mo_map.resolve()), mo_map_dst.name))

    sb_map_dst = stage_dir / "sb_group_map.tsv"
    shutil.copy2(args.sb_group_map, sb_map_dst)
    manifest_rows.append(("map", str(args.sb_group_map.resolve()), sb_map_dst.name))

    rna_solo_dirs = [path.resolve() for path in args.rna_solo_dir]
    rna_filtered_bams = [path.resolve() for path in args.rna_filtered_bam]
    dna_nodup_bams = [path.resolve() for path in args.dna_nodup_bam]

    for solo_dir in rna_solo_dirs:
        dst = stage_dir / solo_dir.name
        copy_directory(solo_dir, dst)
        manifest_rows.append(("rna_solo_dir", str(solo_dir), dst.name))

    for bam in rna_filtered_bams:
        dst = stage_dir / bam.name
        link_or_copy_file(bam, dst)
        manifest_rows.append(("rna_filtered_bam", str(bam), dst.name))

    for bam in dna_nodup_bams:
        dst = stage_dir / bam.name
        link_or_copy_file(bam, dst)
        manifest_rows.append(("dna_nodup_bam", str(bam), dst.name))

        bai = bam.with_name(f"{bam.name}.bai")
        if bai.exists():
            bai_dst = stage_dir / bai.name
            link_or_copy_file(bai, bai_dst)
            manifest_rows.append(("dna_nodup_bai", str(bai), bai_dst.name))

    rna_samples = infer_rna_samples(rna_solo_dirs, rna_filtered_bams)
    dna_group_keys = infer_dna_group_keys(dna_nodup_bams)
    exact_pairs = write_pairs(stage_dir / "pairs.tsv", dna_group_keys, rna_samples)

    write_manifest(stage_dir / "stage_manifest.tsv", manifest_rows)
    write_readiness_note(stage_dir / "sc_process_readiness.txt", args.species, args.genome, exact_pairs)


if __name__ == "__main__":
    main()
