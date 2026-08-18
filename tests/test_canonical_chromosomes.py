import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RESOLVER_PATH = REPO_ROOT / "bin/resolve_canonical_chromosomes.py"
SPEC = importlib.util.spec_from_file_location("canonical_resolver", RESOLVER_PATH)
RESOLVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RESOLVER)


class CanonicalChromosomeTests(unittest.TestCase):
    def test_ucsc_contract_keeps_exact_primary_names_in_reference_order(self):
        entries = [
            ("chr1", 100),
            ("chr1_KI270706v1_random", 20),
            ("chrX", 90),
            ("chrY", 80),
            ("chrM", 70),
            ("chrUn_KI270442v1", 10),
            ("chrEBV", 30),
        ]

        style, canonical = RESOLVER.resolve_canonical_entries(entries, "synthetic UCSC")

        self.assertEqual(style, "ucsc")
        self.assertEqual(
            canonical,
            [("chr1", 100), ("chrX", 90), ("chrY", 80), ("chrM", 70)],
        )

    def test_ensembl_contract_supports_mt_or_m_without_renaming(self):
        for mitochondrial in ("MT", "M"):
            with self.subTest(mitochondrial=mitochondrial):
                entries = [
                    ("1", 100),
                    ("X", 90),
                    ("Y", 80),
                    (mitochondrial, 70),
                    ("GL000220.1", 10),
                    ("KI270728.1", 20),
                ]

                style, canonical = RESOLVER.resolve_canonical_entries(
                    entries, "synthetic Ensembl"
                )

                self.assertEqual(style, "ensembl")
                self.assertEqual(
                    [name for name, _ in canonical], ["1", "X", "Y", mitochondrial]
                )

    def test_mixed_conventions_fail_clearly(self):
        entries = [("chr1", 100), ("X", 90), ("Y", 80), ("MT", 70)]

        with self.assertRaisesRegex(
            RESOLVER.CanonicalChromosomeError,
            "both UCSC-style .* and Ensembl-style",
        ):
            RESOLVER.resolve_canonical_entries(entries, "mixed reference")

    def test_missing_mitochondrial_anchor_fails(self):
        entries = [("chr1", 100), ("chrX", 90), ("chrY", 80)]

        with self.assertRaisesRegex(
            RESOLVER.CanonicalChromosomeError, "missing: chrM"
        ):
            RESOLVER.resolve_canonical_entries(entries, "missing mitochondria")

    def test_bwa_ann_and_chrom_sizes_must_agree(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            ann = directory / "reference.ann"
            sizes = directory / "reference.chrom.sizes"
            ann.write_text(
                "280 5 11\n"
                "0 1 (null)\n0 100 0\n"
                "0 X (null)\n100 60 0\n"
                "0 Y (null)\n160 50 0\n"
                "0 MT (null)\n210 40 0\n"
                "0 GL000220.1 (null)\n250 30 0\n",
                encoding="utf-8",
            )
            sizes.write_text(
                "1\t100\nX\t60\nY\t50\nMT\t40\nGL000220.1\t30\n",
                encoding="utf-8",
            )

            primary = RESOLVER.resolve_canonical_entries(
                RESOLVER.read_bwa_ann(ann), str(ann)
            )
            secondary = RESOLVER.resolve_canonical_entries(
                RESOLVER.read_chrom_sizes(sizes), str(sizes)
            )
            RESOLVER.verify_matching_contracts(primary, secondary, str(ann), str(sizes))

            mismatched = (secondary[0], [("1", 101), *secondary[1][1:]])
            with self.assertRaisesRegex(
                RESOLVER.CanonicalChromosomeError, "names or lengths disagree"
            ):
                RESOLVER.verify_matching_contracts(
                    primary, mismatched, str(ann), "mismatched sizes"
                )

    def test_rna_coverage_does_not_force_ucsc_prefix(self):
        coverage_script = (
            REPO_ROOT / "scripts/core_runtime/RNA_COVERAGE.sh"
        ).read_text(encoding="utf-8")
        workflow_text = (
            REPO_ROOT / "subworkflows/local/rna_core/main.nf"
        ).read_text(encoding="utf-8")

        self.assertNotIn("--outWigReferencesPrefix", coverage_script)
        self.assertIn("meta.canonical_chrom_sizes", workflow_text)


if __name__ == "__main__":
    unittest.main()
