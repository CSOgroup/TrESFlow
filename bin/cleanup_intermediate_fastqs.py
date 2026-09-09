#!/usr/bin/env python3
"""Safely unlink an explicit set of lifecycle-expired FASTQ paths."""

from __future__ import annotations

import argparse
import os
import stat
import sys
from dataclasses import dataclass
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


@dataclass(frozen=True)
class ValidatedTarget:
    path: Path
    parent: Path
    name: str
    device: int | None
    inode: int | None
    mode: int | None
    size: int
    missing: bool


def is_within(path: Path, root: Path) -> bool:
    try:
        return Path(os.path.commonpath((path, root))) == root
    except ValueError:
        return False


def normalize_target(raw_path: str) -> Path:
    if not raw_path or any(ord(character) < 32 for character in raw_path):
        raise ValueError("FASTQ cleanup paths must be non-empty and contain no control characters")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = Path.cwd() / candidate
    return Path(os.path.abspath(candidate))


def validate_targets(work_dir: Path, raw_paths: list[str], mode: str) -> list[ValidatedTarget]:
    work_root = work_dir.resolve(strict=True)
    if not work_root.is_dir():
        raise ValueError(f"Configured work directory is not a directory: {work_root}")
    if not raw_paths:
        raise ValueError("FASTQ cleanup received no target paths")

    targets: list[ValidatedTarget] = []
    seen: set[Path] = set()
    for raw_path in raw_paths:
        candidate = normalize_target(raw_path)
        if candidate in seen:
            raise ValueError(f"Duplicate FASTQ cleanup target: {candidate}")
        seen.add(candidate)

        if not candidate.name.lower().endswith(FASTQ_SUFFIXES):
            raise ValueError(f"Refusing to unlink a non-FASTQ path: {candidate}")

        parent = candidate.parent.resolve(strict=True)
        if not is_within(parent, work_root):
            raise ValueError(
                f"Refusing to unlink a path outside the controlled work directory: {candidate}"
            )

        try:
            metadata = candidate.lstat()
        except FileNotFoundError:
            targets.append(
                ValidatedTarget(candidate, parent, candidate.name, None, None, None, 0, True)
            )
            continue

        if stat.S_ISLNK(metadata.st_mode):
            if mode == "producer-output":
                raise ValueError(f"Refusing producer-output cleanup of a symlink: {candidate}")
            resolved_target = candidate.resolve(strict=True)
            if not is_within(resolved_target, work_root):
                raise ValueError(
                    f"Refusing to unlink a staged symlink whose target is outside the work directory: "
                    f"{candidate} -> {resolved_target}"
                )
            if not resolved_target.is_file():
                raise ValueError(f"Staged FASTQ symlink target is not a regular file: {candidate}")
        elif not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Refusing to unlink a non-regular FASTQ path: {candidate}")

        targets.append(
            ValidatedTarget(
                candidate,
                parent,
                candidate.name,
                metadata.st_dev,
                metadata.st_ino,
                metadata.st_mode,
                metadata.st_size,
                False,
            )
        )

    return targets


def unlink_targets(targets: list[ValidatedTarget]) -> list[tuple[Path, str, int]]:
    results: list[tuple[Path, str, int]] = []
    for target in targets:
        if target.missing:
            results.append((target.path, "already_missing", 0))
            continue

        parent_fd = os.open(target.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            try:
                current = os.stat(target.name, dir_fd=parent_fd, follow_symlinks=False)
            except FileNotFoundError:
                results.append((target.path, "already_missing", 0))
                continue

            if (
                current.st_dev != target.device
                or current.st_ino != target.inode
                or current.st_mode != target.mode
            ):
                raise RuntimeError(f"FASTQ cleanup target changed after validation: {target.path}")
            os.unlink(target.name, dir_fd=parent_fd)
            results.append((target.path, "removed", target.size))
        finally:
            os.close(parent_fd)
    return results


def write_report(report: Path | None, results: list[tuple[Path, str, int]]) -> None:
    if report is None:
        return
    with report.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("path\tstatus\tbytes\n")
        for path, status, size in results:
            handle.write(f"{path}\t{status}\t{size}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument(
        "--mode",
        required=True,
        choices=("producer-output", "staged-input"),
    )
    parser.add_argument("--path", action="append", required=True, dest="paths")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        targets = validate_targets(args.work_dir, args.paths, args.mode)
        results = unlink_targets(targets)
        write_report(args.report, results)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FASTQ cleanup refused or failed: {error}", file=sys.stderr)
        return 1

    if not args.quiet:
        removed = sum(status == "removed" for _, status, _ in results)
        reclaimed = sum(size for _, status, size in results if status == "removed")
        print(f"FASTQ cleanup removed {removed} path(s), {reclaimed} byte(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
