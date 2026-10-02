#!/usr/bin/env python3
"""Validate canonical STAR String identities and remove repeated SAM cell tags.

STAR copies input tags and emits CB alongside corrected UB, including a second
CB when CB was already in the input. A rejected UMI produces CB:Z:-, which is
a placeholder rather than another cell identity. Counting has already used
full XI. This stream preserves alignment order, every alignment class and
all molecular tags, including STAR's UB placeholders.
"""

import sys


def normalize_record(line):
    if line.startswith("@"):
        return line
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        raise ValueError("Malformed STAR SAM record")
    cell = {}
    tags = []
    for token in fields[11:]:
        key = token.split(":", 1)[0]
        if key in ("CB", "XI"):
            if not token.startswith(key + ":Z:"):
                raise ValueError(f"Invalid cell attribute for {fields[0]}: {token}")
            value = token[5:]
            if key == "CB" and value == "-":
                # STAR appends this for rejected UMIs. Require the canonical
                # input CB below and still compare every non-placeholder tag.
                continue
            if key in cell and cell[key] != value:
                raise ValueError(f"Conflicting {key} tags for {fields[0]}")
            cell[key] = value
        else:
            if key in ("CR", "UR") and any(c not in "ACGTN" for c in token[5:]):
                raise ValueError(f"Non-molecular {key} sequence for {fields[0]}")
            tags.append(token)
    if not cell.get("XI") or cell.get("CB") != cell["XI"]:
        raise ValueError(f"STAR CB/XI identity mismatch for {fields[0]}")
    parts = cell["XI"].rsplit("_", 2)
    if len(parts) != 3 or not parts[1].isascii() or not parts[1].isdigit() or int(parts[1]) <= 0:
        raise ValueError(f"Invalid full-cell-v1 identity for {fields[0]}")
    if any(c not in "ACGT" for c in parts[2]) or not parts[2]:
        raise ValueError(f"Invalid ligations for {fields[0]}")
    return "\t".join(fields[:11] + tags + ["CB:Z:" + cell["XI"], "XI:Z:" + cell["XI"]]) + "\n"


def main():
    for line in sys.stdin:
        sys.stdout.write(normalize_record(line))


if __name__ == "__main__":
    main()
