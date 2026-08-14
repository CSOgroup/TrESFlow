#!/usr/bin/env python3

import gzip
import os
import re
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
CANONICAL_CELL_TAG = "XI"
AVITI_QNAME_PATTERN = re.compile(
    r"^[^:]+:[^:]+:[^:]+:([0-9]+):([0-9]+):([0-9]+):([0-9]+):[^:]+$"
)


def resolve_temp_root() -> Path:
    configured = os.environ.get("TMPDIR")
    if not configured:
        raise RuntimeError("TMPDIR is not set. Configure runtime.tmpdir in the samplesheet.")
    root = Path(configured).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    return root.resolve()


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def timestamp() -> str:
    return datetime.now().isoformat(timespec="seconds")


def file_size(path: Path):
    try:
        return path.stat().st_size
    except FileNotFoundError:
        return "missing"


def log_event(message: str, *paths: Path, elapsed: float = None):
    details = [f"{path} size={file_size(path)}" for path in paths]
    if elapsed is not None:
        details.append(f"elapsed={elapsed:.2f}s")
    suffix = f" | {'; '.join(details)}" if details else ""
    sys.stderr.write(f"[{timestamp()}] {message}{suffix}\n")
    sys.stderr.flush()


def fastq_compression_state(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".fastq.gz", ".fq.gz")):
        return "gzipped"
    if name.endswith((".fastq", ".fq")):
        return "uncompressed"
    raise RuntimeError(
        f"Unrecognized FASTQ extension for {path}; expected .fastq.gz, .fq.gz, .fastq, or .fq"
    )


def strict_move_fastq(source: Path, destination: Path):
    source_state = fastq_compression_state(source)
    destination_state = fastq_compression_state(destination)
    if source_state != destination_state:
        raise RuntimeError(
            f"FASTQ compression state mismatch: source={source} ({source_state}) "
            f"destination={destination} ({destination_state}). "
            "Python wrappers must not recompress production FASTQs."
        )

    start = time.monotonic()
    log_event("Starting strict FASTQ move", source, destination)
    shutil.move(source, destination)
    log_event("Finished strict FASTQ move", destination, elapsed=time.monotonic() - start)


def fastq_iter(path: Path):
    with open_maybe_gzip(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def parse_header(header: str):
    if not header.startswith("@"):
        raise ValueError(f"FASTQ header does not start with '@': {header}")
    body = header[1:]
    parts = body.split(" ", 1)
    if len(parts) == 1:
        return parts[0], ""
    return parts[0], parts[1]


def load_whitelist(path: Path):
    with open(path, "rt", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def resolve_codon_bin() -> str:
    configured = os.environ.get("CODON_BIN")
    if configured:
        codon_bin = Path(configured)
        if not codon_bin.exists() or not os.access(codon_bin, os.X_OK):
            raise RuntimeError(f"Configured CODON_BIN is missing or not executable: {codon_bin}")
        return str(codon_bin)

    resolved = shutil.which("codon")
    if resolved is None:
        raise RuntimeError("codon executable not found in PATH")
    return resolved


def find_existing_output(base_dir: Path, candidate_names, label: str) -> Path:
    for candidate_name in candidate_names:
        candidate = base_dir / candidate_name
        if candidate.exists():
            return candidate
    joined = ", ".join(str(base_dir / name) for name in candidate_names)
    raise FileNotFoundError(f"Expected {label} in one of: {joined}")


def stem_without_fastq_suffix(name: str) -> str:
    for suffix in FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def tagged_fastq_candidates(input_name: str, tag: str):
    stem = stem_without_fastq_suffix(input_name)
    return [
        f"{stem}_{tag}.fastq.gz",
        f"{stem}_{tag}.fq.gz",
        f"{stem}_{tag}.fastq",
        f"{stem}_{tag}.fq",
    ]


def normalize_split_fastq_name(name: str) -> str:
    if name.endswith(".fastq.gz"):
        return name
    if name.endswith(".fq.gz"):
        return name[: -len(".fq.gz")] + ".fastq.gz"
    if name.endswith(".fastq"):
        return name
    if name.endswith(".fq"):
        return name[: -len(".fq")] + ".fastq"
    raise RuntimeError(f"Unrecognized split FASTQ extension: {name}")


def move_split_output(source: Path, output_dir: Path) -> Path:
    if source.name.endswith(FASTQ_SUFFIXES):
        destination = output_dir / normalize_split_fastq_name(source.name)
    else:
        destination = output_dir / source.name

    if destination.exists():
        raise RuntimeError(f"Refusing to overwrite existing split output: {destination}")

    start = time.monotonic()
    log_event("Starting split output move", source, destination)
    shutil.move(source, destination)
    log_event("Finished split output move", destination, elapsed=time.monotonic() - start)
    return destination


def percent(count: int, total: int) -> str:
    if total == 0:
        return "0.0%"
    return f"{(count / total) * 100.0}%"


def normalize_sb_drop_first(sb: str):
    if len(sb) < 2:
        raise ValueError(f"SB tag length < 2: {sb}")
    return sb[1:]


def find_tag_value(comment: str, tag_name: str):
    for token in comment.replace("\t", " ").split():
        if token.startswith(f"{tag_name}:"):
            return token.rsplit(":", 1)[-1]
    return ""


def canonical_cell_id(sample: str, group_name: str, cell_barcode: str) -> str:
    return f"{sample}_{group_name}_{cell_barcode}"


def parse_aviti_qname(read_name: str):
    """Return AVITI lane, tile, x, and y as integers from a FASTQ/BAM read name."""
    match = AVITI_QNAME_PATTERN.fullmatch(read_name)
    if match is None:
        raise ValueError(
            "DNA read name does not match the expected AVITI format "
            f"instrument:run:flowcell:lane:tile:x:y:UMI: {read_name}"
        )
    return tuple(int(value) for value in match.groups())


def aviti_lane_read_group(lane: int) -> str:
    """Build the small SAM-safe read-group ID used for an AVITI lane."""
    return f"L{lane}"


def cell_barcode_without_sb(cb: str, sb: str, sample: str, group_name: str) -> str:
    candidate_prefixes = [sb]
    if len(sb) > 1:
        candidate_prefixes.append(sb[1:])

    for prefix in candidate_prefixes:
        if cb.startswith(prefix):
            cell_barcode = cb[len(prefix) :]
            if cell_barcode:
                return cell_barcode

    raise ValueError(
        f"Cannot derive canonical cell barcode for sample {sample} group {group_name}: "
        f"CB tag '{cb}' does not start with SB tag '{sb}'"
    )


def canonicalize_fastq_comment(
    sample: str,
    group_name: str,
    comment: str,
    read_group: str = None,
) -> str:
    cb = find_tag_value(comment, "CB")
    sb = find_tag_value(comment, "SB")
    if not cb or not sb:
        raise ValueError(
            f"Missing CB or SB tag while canonicalizing cell ID for sample {sample} group {group_name}"
        )

    technical_cell = cell_barcode_without_sb(cb, sb, sample, group_name)
    output_read_group = read_group if read_group is not None else technical_cell
    canonical = canonical_cell_id(sample, group_name, technical_cell)
    tokens = []
    has_rg = False
    has_canonical = False
    for token in comment.replace("\t", " ").split():
        if token.startswith("CB:"):
            tokens.append(f"CB:Z:{technical_cell}")
        elif token.startswith("RG:"):
            tokens.append(f"RG:Z:{output_read_group}")
            has_rg = True
        elif token.startswith(f"{CANONICAL_CELL_TAG}:"):
            tokens.append(f"{CANONICAL_CELL_TAG}:Z:{canonical}")
            has_canonical = True
        else:
            tokens.append(token)
    if not has_rg:
        tokens.append(f"RG:Z:{output_read_group}")
    if not has_canonical:
        tokens.append(f"{CANONICAL_CELL_TAG}:Z:{canonical}")
    return "\t".join(tokens)


def canonicalize_dna_fastq_comment(
    sample: str,
    group_name: str,
    read_name: str,
    comment: str,
) -> str:
    """Canonicalize DNA cell tags while making RG identify the AVITI lane."""
    lane, _, _, _ = parse_aviti_qname(read_name)
    cb = find_tag_value(comment, "CB")
    sb = find_tag_value(comment, "SB")
    if not cb or not sb:
        raise ValueError(
            f"Missing CB or SB tag while canonicalizing cell ID for sample {sample} group {group_name}"
        )
    technical_cell = cell_barcode_without_sb(cb, sb, sample, group_name)
    read_group = aviti_lane_read_group(lane)
    return canonicalize_fastq_comment(sample, group_name, comment, read_group=read_group)


def load_sb_group_map(path: Path, sample: str):
    sb_to_group = {}
    group_names = []
    group_seen = set()

    with open(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            row_sample, group_name, sb_bc = parts[0], parts[1], parts[2]
            if row_sample != sample:
                continue

            if group_name not in group_seen:
                group_names.append(group_name)
                group_seen.add(group_name)

            if sb_bc in sb_to_group and sb_to_group[sb_bc] != group_name:
                raise ValueError(f"SB group conflict for sample {sample} SB {sb_bc}")
            sb_to_group[sb_bc] = group_name

    if not sb_to_group:
        raise ValueError(f"No SB group mapping found for sample {sample} in {path}")

    return sb_to_group, group_names


def resolve_group(sample: str, sb_raw: str, sb_to_group):
    if sb_raw in sb_to_group:
        return sb_to_group[sb_raw]

    key = normalize_sb_drop_first(sb_raw)
    if key in sb_to_group:
        return sb_to_group[key]

    raise ValueError(f"SB not found in SB group map for sample {sample}: raw={sb_raw} key={key}")


def write_fastq_record(handle, name: str, comment: str, seq: str, qual: str):
    handle.write("@")
    handle.write(name)
    if comment:
        handle.write(" ")
        handle.write(comment)
    handle.write("\n")
    handle.write(seq)
    handle.write("\n+\n")
    handle.write(qual)
    handle.write("\n")


def write_rg_header(path: Path, sample: str, library_name: str, read_groups):
    with open(path, "wt", encoding="utf-8") as handle:
        for read_group in sorted(read_groups):
            # All lane RGs share LB so Picard still groups genomic duplicates
            # across lanes; CB provides cell identity and RG only separates lanes.
            handle.write(
                f"@RG\tID:{read_group}\tSM:{sample}\tLB:{library_name}"
                "\tPL:ELEMENT\tPM:AVITI_500MIO\n"
            )
