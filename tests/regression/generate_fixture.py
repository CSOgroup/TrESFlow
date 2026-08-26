#!/usr/bin/env python3
"""Build deterministic, test-only TrESFlow inputs and miniature references.

The generated assembly is artificial and contains no biological source data.
STAR and BWA-MEM2 indices are deliberately written outside the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


SCHEMA_VERSION = 1
GENERATOR_VERSION = "1.0.0"
ASSEMBLY_ID = "tresflow-phase0-synthetic-v1"
ASSEMBLY_LICENSE = "CC0-1.0"
EXPECTED_STAR_VERSION = "2.7.11b"
EXPECTED_BWA_MEM2_VERSION = "2.2.1"
CONTIG_LENGTHS = {"chr1": 50_000, "chrX": 4_000, "chrY": 4_000, "chrM": 4_000}
STAR_ARGUMENTS = [
    "--runMode", "genomeGenerate",
    "--runThreadN", "2",
    "--sjdbOverhang", "99",
    "--genomeSAindexNbases", "6",
    "--genomeChrBinNbits", "12",
]

RNA_I1_TEMPLATE = (
    "AAAAAAAAAAAAAAAACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGCATGCA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGATCGATCAAAAAAAAAAAAAAAAAAAAA"
)
DNA_SINGLE_I1_TEMPLATE = "A" * 15 + "ACGTACGT" + "C" * 30 + "TGCATGCA" + "C" * 30 + "GATCGATC"
DNA_SINGLE_I2_TEMPLATE = "A" * 14 + "TACG" + "TACGTAAA" + "C" * 10
DNA_DUAL_I1_TEMPLATE = "AAA" + "AGGCTATA" + "C" * 30 + "ACGTACGT" + "C" * 30 + "TGCATGCA" + "C" * 30 + "GATCGATC"
DUAL_TAG_ARTIFACT = "AAGTATGCAGCGCGCTCAAGCAC"


def deterministic_sequence(label: str, length: int) -> str:
    """Return a stable pseudo-random DNA sequence without RNG-version coupling."""
    sequence: list[str] = []
    counter = 0
    while len(sequence) < length:
        digest = hashlib.sha256(f"{ASSEMBLY_ID}:{label}:{counter}".encode()).digest()
        sequence.extend("ACGT"[byte & 3] for byte in digest)
        counter += 1
    return "".join(sequence[:length])


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def encode_umi(number: int) -> str:
    bases = "ACGT"
    chars = []
    for _ in range(10):
        chars.append(bases[number & 3])
        number >>= 2
    return "".join(reversed(chars))


def fasta_text(contigs: dict[str, str]) -> str:
    lines: list[str] = []
    for name, sequence in contigs.items():
        lines.append(f">{name}")
        lines.extend(sequence[offset : offset + 80] for offset in range(0, len(sequence), 80))
    return "\n".join(lines) + "\n"


def fastq_text(records: list[tuple[str, str]]) -> str:
    lines: list[str] = []
    for name, sequence in records:
        lines.extend((f"@{name}", sequence, "+", "I" * len(sequence)))
    return "\n".join(lines) + "\n"


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_version(command: Path, args: list[str]) -> str:
    result = subprocess.run(
        [str(command), *args], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        raise SystemExit(f"Version command produced no output: {command} {' '.join(args)}")
    # BWA-MEM2 wrappers can emit SIMD launcher messages before the version.
    return lines[-1]


def build_fastqs(root: Path, contigs: dict[str, str], count: int) -> dict[str, dict[str, Path]]:
    reads_root = root / "reads"

    rna_r1_positive = contigs["chr1"][2_000:2_100]
    rna_r2_positive = reverse_complement(contigs["chr1"][2_300:2_386])
    rna_r1_negative = reverse_complement(contigs["chr1"][5_000:5_100])
    rna_r2_negative = contigs["chr1"][4_700:4_786]
    rna_i1: list[tuple[str, str]] = []
    rna_r1: list[tuple[str, str]] = []
    rna_r2: list[tuple[str, str]] = []
    for index in range(count):
        name = f"phase0_rna_{index:04d}/1"
        mate_name = f"phase0_rna_{index:04d}/2"
        rna_i1.append((name, RNA_I1_TEMPLATE))
        negative_strand = index >= count // 2
        rna_r1.append((name, rna_r1_negative if negative_strand else rna_r1_positive))
        # ACTG reverse-complements to the configured CAGT sample barcode.
        genomic_r2 = rna_r2_negative if negative_strand else rna_r2_positive
        rna_r2.append((mate_name, "ACTG" + encode_umi(index) + genomic_r2))

    dna_r1_sequence = contigs["chr1"][10_000:10_100]
    dna_r2_sequence = reverse_complement(contigs["chr1"][10_300:10_400])

    def dna_records(prefix: str, i1_template: str, i2_template: str, dual: bool) -> dict[str, Path]:
        i1: list[tuple[str, str]] = []
        i2: list[tuple[str, str]] = []
        r1: list[tuple[str, str]] = []
        r2: list[tuple[str, str]] = []
        for index in range(count):
            # Split_ReadsV2 preserves AVITI flow-cell coordinates for optical
            # duplicate handling and requires the native eight-field shape.
            base_name = (
                f"AVTEST:RUN1:FLOW1:1:1101:{1_000 + index}:{2_000 + index}:"
                f"{encode_umi(index)}"
            )
            i1.append((f"{base_name}/1", i1_template))
            i2.append((f"{base_name}/2", i2_template))
            first = dna_r1_sequence
            if dual and index == 0:
                # One pair exercises the exact dual-tag artifact rejection path.
                first = (
                    dna_r1_sequence[:10]
                    + DUAL_TAG_ARTIFACT
                    + dna_r1_sequence[10 + len(DUAL_TAG_ARTIFACT) :]
                )
            r1.append((f"{base_name}/1", first))
            r2.append((f"{base_name}/2", dna_r2_sequence))
        scenario_root = reads_root / prefix
        paths = {
            "i1": scenario_root / "I1.fastq",
            "i2": scenario_root / "I2.fastq",
            "r1": scenario_root / "R1.fastq",
            "r2": scenario_root / "R2.fastq",
        }
        for role, records in (("i1", i1), ("i2", i2), ("r1", r1), ("r2", r2)):
            write_text(paths[role], fastq_text(records))
        return paths

    rna_paths = {
        "i1": reads_root / "rna_only" / "I1.fastq",
        "r1": reads_root / "rna_only" / "R1.fastq",
        "r2": reads_root / "rna_only" / "R2.fastq",
    }
    for role, records in (("i1", rna_i1), ("r1", rna_r1), ("r2", rna_r2)):
        write_text(rna_paths[role], fastq_text(records))

    return {
        "rna_only": rna_paths,
        "dna_single": dna_records(
            "dna_single", DNA_SINGLE_I1_TEMPLATE, DNA_SINGLE_I2_TEMPLATE, False
        ),
        "dna_dual": dna_records(
            "dna_dual", DNA_DUAL_I1_TEMPLATE, DNA_DUAL_I1_TEMPLATE, True
        ),
    }


def samplesheet_text(
    scenario: str,
    paths: dict[str, Path],
    reference_root: Path,
    env_prefix: Path,
    runtime_tmpdir: Path,
) -> str:
    common = f"""library_name: PHASE0_SYNTHETIC

runtime:
  env_prefix: {env_prefix}
  tmpdir: {runtime_tmpdir}

references:
  species: human
  root: {reference_root}
  ligation_barcode_whitelist: {reference_root / 'ligation_barcode_whitelist.txt'}
"""
    if scenario == "rna_only":
        return common + f"""  rna_ref_dir: {reference_root / 'rna/human/star'}

samples:
  phase0_rna:
    groups:
      Normal:
        sb_barcodes: [CAGT]
    rna:
      reads:
        i1: {paths['i1']}
        r1: {paths['r1']}
        r2: {paths['r2']}
"""

    tagmentation = "single" if scenario == "dna_single" else "dual"
    sample_id = f"phase0_{scenario}"
    barcode_key = "sb_barcodes" if scenario == "dna_single" else "dna_sb_barcodes"
    barcode = "CGTA" if scenario == "dna_single" else "AAA"
    i2_line = f"        i2: {paths['i2']}\n" if scenario == "dna_single" else ""
    mark_barcode = "TTTACGTA" if scenario == "dna_single" else "AGGCTATA"
    return common + f"""  dna_ref_dir: {reference_root / 'dna/human/bwa'}
  dna_blacklist_bed: {reference_root / 'dna/human/blacklist.bed'}
  dna_chrom_sizes: {reference_root / 'dna/human/chrom.sizes'}
  dna_effective_genome_size: {sum(CONTIG_LENGTHS.values())}

samples:
  {sample_id}:
    groups:
      Normal:
        {barcode_key}: [{barcode}]
        mark_barcodes:
          H3K27ac: {mark_barcode}
    dna:
      tagmentation: {tagmentation}
      reads:
        i1: {paths['i1']}
{i2_line}        r1: {paths['r1']}
        r2: {paths['r2']}
"""


def build_fixture(args: argparse.Namespace) -> None:
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty fixture directory: {output}")
    output.mkdir(parents=True, exist_ok=True)

    star = args.star.resolve()
    bwa_mem2 = args.bwa_mem2.resolve()
    for executable in (star, bwa_mem2):
        if not executable.is_file():
            raise SystemExit(f"Required executable not found: {executable}")

    star_version = command_version(star, ["--version"])
    bwa_mem2_version = command_version(bwa_mem2, ["version"])
    if star_version != EXPECTED_STAR_VERSION:
        raise SystemExit(
            f"STAR {EXPECTED_STAR_VERSION} is required to reproduce the Phase 0 fixture; "
            f"observed {star_version}"
        )
    if bwa_mem2_version != EXPECTED_BWA_MEM2_VERSION:
        raise SystemExit(
            f"BWA-MEM2 {EXPECTED_BWA_MEM2_VERSION} is required to reproduce the Phase 0 fixture; "
            f"observed {bwa_mem2_version}"
        )

    contigs = {name: deterministic_sequence(name, length) for name, length in CONTIG_LENGTHS.items()}
    reference_root = output / "reference"
    source_root = reference_root / "source"
    fasta = source_root / "genome.fa"
    gtf = source_root / "genes.gtf"
    write_text(fasta, fasta_text(contigs))
    write_text(
        gtf,
        "chr1\tTrESFlow\tgene\t1001\t20000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"SyntheticGene1\";\n"
        "chr1\tTrESFlow\ttranscript\t1001\t20000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; gene_name \"SyntheticGene1\";\n"
        "chr1\tTrESFlow\texon\t1001\t20000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; gene_name \"SyntheticGene1\";\n",
    )
    write_text(reference_root / "ligation_barcode_whitelist.txt", "ACGTACGT\nTGCATGCA\nGATCGATC\n")
    write_text(reference_root / "dna/human/blacklist.bed", "chr1\t30000\t30100\n")
    write_text(
        reference_root / "dna/human/chrom.sizes",
        "".join(f"{name}\t{length}\n" for name, length in CONTIG_LENGTHS.items()),
    )

    star_dir = reference_root / "rna/human/star"
    star_dir.mkdir(parents=True, exist_ok=True)
    star_cmd = [
        str(star),
        *STAR_ARGUMENTS[:4],
        "--genomeDir", str(star_dir),
        "--genomeFastaFiles", str(fasta),
        "--sjdbGTFfile", str(gtf),
        *STAR_ARGUMENTS[4:],
    ]
    subprocess.run(star_cmd, check=True, cwd=output)

    bwa_dir = reference_root / "dna/human/bwa"
    bwa_dir.mkdir(parents=True, exist_ok=True)
    bwa_fasta = bwa_dir / "genome.fa"
    shutil.copyfile(fasta, bwa_fasta)
    bwa_cmd = [str(bwa_mem2), "index", str(bwa_fasta)]
    subprocess.run(bwa_cmd, check=True, cwd=output)

    read_paths = build_fastqs(output, contigs, args.read_pairs)
    sheets_root = output / "samplesheets"
    for scenario, paths in read_paths.items():
        write_text(
            sheets_root / f"{scenario}.yaml",
            samplesheet_text(
                scenario,
                paths,
                reference_root,
                args.env_prefix.resolve(),
                output / "tmp" / scenario,
            ),
        )

    files = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "asset-manifest.json":
            files.append(
                {
                    "path": path.relative_to(output).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "fixture_generator": {
            "path": "tests/regression/generate_fixture.py",
            "version": GENERATOR_VERSION,
            "sha256": sha256(Path(__file__).resolve()),
        },
        "assembly": {
            "id": ASSEMBLY_ID,
            "kind": "deterministic artificial test sequence",
            "license": ASSEMBLY_LICENSE,
            "contigs": CONTIG_LENGTHS,
            "annotation": "one artificial gene/transcript/exon on chr1",
            "chromosome_naming": "UCSC-style test names",
            "effective_genome_size": sum(CONTIG_LENGTHS.values()),
        },
        "read_pairs_per_scenario": args.read_pairs,
        "tools": {
            "STAR": {
                "version": star_version,
                "command": star_cmd,
            },
            "bwa-mem2": {
                "version": bwa_mem2_version,
                "command": bwa_cmd,
            },
        },
        "scenarios": {
            scenario: {
                "samplesheet": f"samplesheets/{scenario}.yaml",
                "modality": "RNA" if scenario == "rna_only" else "DNA",
                "tagmentation": None if scenario == "rna_only" else scenario.removeprefix("dna_"),
                "mock_processes": [],
            }
            for scenario in read_paths
        },
        "files": files,
    }
    write_text(output / "asset-manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--env-prefix", required=True, type=Path)
    parser.add_argument("--star", required=True, type=Path)
    parser.add_argument("--bwa-mem2", required=True, type=Path)
    parser.add_argument("--read-pairs", type=int, default=64)
    args = parser.parse_args()
    if args.read_pairs < 32:
        parser.error("--read-pairs must be at least 32 so STARsolo cell calling is exercised")
    return args


if __name__ == "__main__":
    build_fixture(parse_args())
