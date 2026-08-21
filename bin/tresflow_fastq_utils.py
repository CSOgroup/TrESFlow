#!/usr/bin/env python3

import gzip
import os
import re
import shutil
import sys
import time
from datetime import datetime
from itertools import zip_longest
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
CANONICAL_CELL_TAG = "XI"
AVITI_QNAME_PATTERN = re.compile(
    r"^([^:]+):([^:]+):([^:]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+):([^:]+)$"
)
SAM_SAFE_PHYSICAL_UNIT_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")
CONTROL_CHARACTER_PATTERN = re.compile(r"[\x00-\x1f\x7f]")


def resolve_temp_root() -> Path:
    configured = os.environ.get("TMPDIR")
    if not configured:
        raise RuntimeError("TMPDIR is not set. Configure runtime.tmpdir in the samplesheet.")
    root = Path(configured).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    return root.resolve()


def open_maybe_gzip(path: Path, mode: str):
    if path.name.lower().endswith(".gz"):
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
        record_number = 0
        while True:
            header = handle.readline()
            if header == "":
                break
            record_number += 1
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if seq == "" or plus == "" or qual == "":
                raise ValueError(
                    f"Incomplete FASTQ record {record_number} in {path}: expected four lines"
                )
            header = header.rstrip("\r\n")
            sequence = seq.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = qual.rstrip("\r\n")
            if not header.startswith("@"):
                raise ValueError(
                    f"Malformed FASTQ record {record_number} in {path}: header does not start with '@'"
                )
            if not plus.startswith("+"):
                raise ValueError(
                    f"Malformed FASTQ record {record_number} in {path}: separator does not start with '+'"
                )
            if len(sequence) != len(quality):
                raise ValueError(
                    f"Malformed FASTQ record {record_number} in {path}: sequence and quality lengths differ "
                    f"({len(sequence)} != {len(quality)})"
                )
            yield header, sequence, plus, quality


def read_fastq_manifest(path: Path):
    manifest = path.resolve()
    paths = []
    with open(manifest, "rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            entry = raw_line.rstrip("\r\n")
            if not entry:
                raise ValueError(f"Empty FASTQ manifest entry at {manifest}:{line_number}")
            if CONTROL_CHARACTER_PATTERN.search(entry):
                raise ValueError(f"Control character in FASTQ manifest entry at {manifest}:{line_number}")
            candidate = Path(entry)
            if not candidate.is_absolute():
                candidate = manifest.parent / candidate
            if not candidate.is_file():
                raise ValueError(
                    f"FASTQ manifest entry is not a regular file at {manifest}:{line_number}: {candidate}"
                )
            fastq_compression_state(candidate)
            paths.append(candidate.resolve())
    if not paths:
        raise ValueError(f"FASTQ manifest contains no paths: {manifest}")
    if len(set(paths)) != len(paths):
        raise ValueError(f"FASTQ manifest contains duplicate paths: {manifest}")
    return paths


def resolve_fastq_paths(single_path: Path = None, manifest_path: Path = None):
    if (single_path is None) == (manifest_path is None):
        raise ValueError("Specify exactly one FASTQ path or FASTQ manifest")
    if manifest_path is not None:
        return read_fastq_manifest(manifest_path)
    if not single_path.is_file():
        raise ValueError(f"FASTQ input is not a regular file: {single_path}")
    fastq_compression_state(single_path)
    return [single_path.resolve()]


def fastq_input_spec(single_path: Path = None, manifest_path: Path = None) -> str:
    if manifest_path is not None:
        return f"manifest:{manifest_path.resolve()}"
    return str(single_path.resolve())


def normalize_qname(name: str) -> str:
    if name.endswith(("/1", "/2")):
        return name[:-2]
    return name


def read_read_set_counts(path: Path):
    counts = []
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n")
        if header != "read_set_index\trecord_count":
            raise ValueError(f"Invalid technical read-set count header in {path}: {header!r}")
        for line_number, raw_line in enumerate(handle, start=2):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0].isdigit() or not fields[1].isdigit():
                raise ValueError(f"Invalid technical read-set count row at {path}:{line_number}")
            expected_index = len(counts) + 1
            if int(fields[0]) != expected_index:
                raise ValueError(
                    f"Non-sequential technical read-set index at {path}:{line_number}: "
                    f"expected {expected_index}, found {fields[0]}"
                )
            counts.append(int(fields[1]))
    if not counts:
        raise ValueError(f"Technical read-set count file contains no rows: {path}")
    return counts


def write_read_set_counts(path: Path, counts):
    lines = ["read_set_index\trecord_count"]
    lines.extend(f"{index}\t{count}" for index, count in enumerate(counts, start=1))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _validate_synchronized_qnames(roles, records, read_set_index, record_number):
    qnames = {}
    for role, record in zip(roles, records):
        name, _ = parse_header(record[0])
        qnames[role] = normalize_qname(name)
    if len(set(qnames.values())) != 1:
        detail = ", ".join(f"{role}={name}" for role, name in qnames.items())
        raise ValueError(
            f"Synchronized FASTQ QNAME mismatch in technical read set {read_set_index + 1} "
            f"at record {record_number}: {detail}"
        )


def synchronized_fastq_iter(role_paths, expected_read_set_counts=None, observed_counts=None):
    """Yield complete synchronized records across ordered virtual FASTQ streams."""
    roles = list(role_paths)
    path_counts = {role: len(role_paths[role]) for role in roles}
    if expected_read_set_counts is None and len(set(path_counts.values())) != 1:
        detail = ", ".join(f"{role}={path_counts[role]}" for role in roles)
        raise ValueError(f"Synchronized FASTQ manifests have conflicting read-set counts: {detail}")

    if expected_read_set_counts is None:
        read_set_count = path_counts[roles[0]]
        combined_iterators = {}
    else:
        read_set_count = len(expected_read_set_counts)
        invalid = [role for role in roles if path_counts[role] not in (1, read_set_count)]
        if invalid:
            detail = ", ".join(f"{role}={path_counts[role]}" for role in roles)
            raise ValueError(
                f"FASTQ manifests must contain either one combined path or {read_set_count} "
                f"technical read-set paths when boundaries are supplied: {detail}"
            )
        combined_iterators = {
            role: fastq_iter(role_paths[role][0]) for role in roles if path_counts[role] == 1
        }

    for read_set_index in range(read_set_count):
        iterators = [
            combined_iterators[role]
            if role in combined_iterators
            else fastq_iter(role_paths[role][read_set_index])
            for role in roles
        ]
        record_count = 0
        if expected_read_set_counts is None:
            record_source = enumerate(zip_longest(*iterators), start=1)
        else:
            expected_count = expected_read_set_counts[read_set_index]
            record_source = (
                (record_number, tuple(next(iterator, None) for iterator in iterators))
                for record_number in range(1, expected_count + 1)
            )

        for record_number, records in record_source:
            exhausted = [role for role, record in zip(roles, records) if record is None]
            if exhausted:
                if expected_read_set_counts is None and len(exhausted) == len(roles):
                    break
                remaining = [role for role in roles if role not in exhausted]
                raise ValueError(
                    f"Synchronized FASTQ streams have unequal EOF in technical read set "
                    f"{read_set_index + 1} at record {record_number}: "
                    f"exhausted={','.join(exhausted)}; remaining={','.join(remaining)}"
                )
            record_count += 1
            _validate_synchronized_qnames(roles, records, read_set_index, record_number)
            yield records

        if expected_read_set_counts is not None:
            for role, iterator in zip(roles, iterators):
                if role not in combined_iterators and next(iterator, None) is not None:
                    raise ValueError(
                        f"Synchronized FASTQ stream {role} has more than "
                        f"{expected_read_set_counts[read_set_index]} records in technical read set "
                        f"{read_set_index + 1}"
                    )

        if observed_counts is not None:
            observed_counts.append(record_count)

    for role, iterator in combined_iterators.items():
        if next(iterator, None) is not None:
            expected_total = sum(expected_read_set_counts)
            raise ValueError(
                f"Combined synchronized FASTQ stream {role} has more than {expected_total} records "
                "defined by the technical read-set boundaries"
            )


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
            fields = token.split(":", 2)
            return fields[2] if len(fields) == 3 else ""
    return ""


def canonical_cell_id(sample: str, group_name: str, cell_barcode: str) -> str:
    return f"{sample}_{group_name}_{cell_barcode}"


def parse_aviti_qname(read_name: str):
    """Return AVITI physical-unit fields plus tile/x/y coordinates."""
    match = AVITI_QNAME_PATTERN.fullmatch(read_name)
    if match is None:
        raise ValueError(
            "DNA read name does not match the expected AVITI format "
            f"instrument:run:flowcell:lane:tile:x:y:UMI: {read_name}"
        )
    instrument, run, flowcell, lane, tile, x, y, _umi = match.groups()
    for label, value in (
        ("instrument", instrument),
        ("run", run),
        ("flowcell", flowcell),
    ):
        if SAM_SAFE_PHYSICAL_UNIT_PATTERN.fullmatch(value) is None:
            raise ValueError(
                f"Unsupported character in AVITI {label} for SAM read-group ID: {value!r} ({read_name})"
            )
    return instrument, run, flowcell, int(lane), int(tile), int(x), int(y)


def aviti_physical_unit(instrument: str, run: str, flowcell: str, lane: int) -> str:
    """Build a deterministic SAM-safe AVITI physical-unit/read-group ID."""
    return f"{instrument}:{run}:{flowcell}:L{lane}"


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
    """Canonicalize DNA cell tags while making RG identify its physical unit."""
    instrument, run, flowcell, lane, _, _, _ = parse_aviti_qname(read_name)
    cb = find_tag_value(comment, "CB")
    sb = find_tag_value(comment, "SB")
    if not cb or not sb:
        raise ValueError(
            f"Missing CB or SB tag while canonicalizing cell ID for sample {sample} group {group_name}"
        )
    technical_cell = cell_barcode_without_sb(cb, sb, sample, group_name)
    read_group = aviti_physical_unit(instrument, run, flowcell, lane)
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
            # All physical-unit RGs share LB so Picard still groups genomic
            # duplicates across units; CB provides cell identity.
            handle.write(
                f"@RG\tID:{read_group}\tSM:{sample}\tLB:{library_name}\tPU:{read_group}"
                "\tPL:ELEMENT\tPM:AVITI_500MIO\n"
            )
