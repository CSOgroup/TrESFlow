import csv
import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location(
    "barcode_gate_metrics", REPO_ROOT / "bin/write_barcode_gate_metrics.py"
)
METRICS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(METRICS)


def write_fastqs(root, names):
    for mate in (1, 2):
        path = root / f"R{mate}.fastq"
        with path.open("w", encoding="utf-8") as handle:
            for name in names:
                handle.write(f"@{name}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")


def write_records(path, rows):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for name, l1, l2, l3, sb, mo in rows:
            fields = [name, "CB:Z:recorded", f"L1:Z:{l1}", f"L2:Z:{l2}", f"L3:Z:{l3}", f"SB:Z:{sb}"]
            if mo is not None:
                fields.append(f"MO:Z:{mo}")
            handle.write("\t".join(fields) + "\n")


def read_tsv(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class BarcodeGateMetricTests(unittest.TestCase):
    def run_metrics(self, modality, retained_names, records):
        root_context = tempfile.TemporaryDirectory()
        root = Path(root_context.name)
        write_fastqs(root, retained_names)
        tag_records = root / "records.tsv.gz"
        write_records(tag_records, records)
        sb_map = root / "sb.tsv"
        sb_map.write_text("sample\tgroup1\tAAA\nsample\tgroup2\tCCC\n", encoding="utf-8")
        mo_map = root / "mo.tsv"
        mo_map.write_text(
            "sample\tgroup1\tH3K27ac\tMARKA\n"
            "sample\tgroup1\tH3K27me3\tMARKB\n"
            "sample\tgroup2\tH3K9me3\tMARKC\n",
            encoding="utf-8",
        )
        argv = [
            "--sample", "sample", "--modality", modality,
            "--r1", str(root / "R1.fastq"), "--r2", str(root / "R2.fastq"),
            "--tag-records", str(tag_records), "--sb-group-map", str(sb_map),
            "--output-gates", str(root / "gates.tsv"),
            "--output-composition", str(root / "composition.tsv"),
        ]
        if modality == "dna":
            argv.extend(["--mo-map", str(mo_map)])
        METRICS.main(argv)
        return root_context, read_tsv(root / "gates.tsv"), read_tsv(root / "composition.tsv")

    def test_rna_uses_exact_same_pair_intersections_and_exhaustive_composition(self):
        records = [
            ("raw_removed", "A", "B", "C", "AAA", None),
            ("q1", "A", "B", "C", "AAA", None),
            ("q2", "A", "B", "NoMatch", "CCC", None),
            ("q3", "A", "B", "C", "NoMatch", None),
            ("q4", "A", "B", "C", "CCC", None),
        ]
        context, gates, composition = self.run_metrics("rna", ["q1", "q2", "q3", "q4"], records)
        try:
            values = {row["metric"]: int(row["pairs"]) for row in gates}
            self.assertEqual(
                values,
                {
                    "split_input_pairs": 4,
                    "ligation_barcode_accepted_pairs": 3,
                    "sample_barcode_accepted_pairs": 2,
                },
            )
            self.assertGreaterEqual(values["split_input_pairs"], values["ligation_barcode_accepted_pairs"])
            self.assertGreaterEqual(values["ligation_barcode_accepted_pairs"], values["sample_barcode_accepted_pairs"])
            counts = {row["category"]: int(row["count"]) for row in composition}
            self.assertEqual(counts, {"group1": 1, "group2": 1, "NoMatch": 1})
            self.assertEqual(sum(counts.values()), 3)
        finally:
            context.cleanup()

    def test_dna_gates_are_nested_and_mark_categories_reconcile_per_run(self):
        records = [
            ("q1", "A", "B", "C", "AAA", "MARKA"),
            ("q2", "A", "B", "C", "AAA", "NoMatch"),
            ("q3", "A", "B", "C", "CCC", "MARKC"),
            ("q4", "A", "B", "NoMatch", "AAA", "MARKB"),
            ("q5", "A", "B", "C", "NoMatch", "MARKA"),
        ]
        context, gates, composition = self.run_metrics("dna", [f"q{i}" for i in range(1, 6)], records)
        try:
            values = {row["metric"]: int(row["pairs"]) for row in gates}
            self.assertEqual(list(values.values()), [5, 4, 3, 2])
            sample_rows = [row for row in composition if row["barcode_type"] == "sample_barcode"]
            self.assertEqual(sum(int(row["count"]) for row in sample_rows), 4)
            group1 = [row for row in composition if row["barcode_type"] == "dna_mark" and row["group"] == "group1"]
            group2 = [row for row in composition if row["barcode_type"] == "dna_mark" and row["group"] == "group2"]
            self.assertEqual({row["category"]: int(row["count"]) for row in group1}, {"H3K27ac": 1, "H3K27me3": 0, "NoMatch": 1})
            self.assertEqual({row["category"]: int(row["count"]) for row in group2}, {"H3K9me3": 1, "NoMatch": 0})
            self.assertEqual(sum(int(row["count"]) for row in group1), 2)
            self.assertEqual(sum(int(row["count"]) for row in group2), 1)
        finally:
            context.cleanup()


if __name__ == "__main__":
    unittest.main()
