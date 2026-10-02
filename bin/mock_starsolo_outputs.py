#!/usr/bin/env python3
"""Small test-profile STARsolo outputs retaining the real identity contract."""

import sys
from collections import defaultdict
from pathlib import Path


def main():
    usam, output = map(Path, sys.argv[1:])
    cells = defaultdict(set)
    reads = defaultdict(int)
    for line in usam.read_text().splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        if not int(fields[1]) & 0x40:
            continue
        tags = dict(token.split(":", 2)[::2] for token in fields[11:])
        if tags["CB"] != tags["XI"]:
            raise ValueError("Mock input CB/XI mismatch")
        cells[tags["XI"]].add(tags["UM"])
        reads[tags["XI"]] += 1
    barcodes = sorted(cells)
    matrix = "%%MatrixMarket matrix coordinate integer general\n%\n"
    matrix += f"1 {len(barcodes)} {len(barcodes)}\n"
    matrix += "".join(f"1 {i} {len(cells[cb])}\n" for i, cb in enumerate(barcodes, 1))
    for layer in ("raw", "filtered"):
        directory = output / layer
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "barcodes.tsv").write_text("".join(cb + "\n" for cb in barcodes))
        (directory / "features.tsv").write_text("mock_feature\tmock_feature\tGene Expression\n")
        (directory / "matrix.mtx").write_text(matrix)
    (output / "CellReads.stats").write_text("CB\tcbMatch\tnUMIunique\n" + "".join(
        f"{cb}\t{reads[cb]}\t{len(cells[cb])}\n" for cb in barcodes))
    total_reads = sum(reads.values())
    total_umis = sum(len(values) for values in cells.values())
    saturation = 1 - total_umis / total_reads if total_reads else 0
    (output / "Summary.csv").write_text(
        f"Number of Reads,{total_reads}\n"
        f"Sequencing Saturation,{saturation}\n"
        "Reads Mapped to Genome: Unique+Multiple,1.00\n"
        "Reads Mapped to Genome: Unique,1.00\n"
        "Reads Mapped to GeneFull: Unique+Multiple GeneFull,1.00\n"
        "Reads Mapped to GeneFull: Unique GeneFull,1.00\n"
        f"Estimated Number of Cells,{len(barcodes)}\n"
        f"UMIs in Cells,{total_umis}\n"
    )


if __name__ == "__main__":
    main()
