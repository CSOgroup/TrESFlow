#!/usr/bin/env python3

import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


fixture = load_module("phase0_fixture", "tests/regression/generate_fixture.py")
normalizer = load_module("phase0_normalizer", "tests/regression/normalize_outputs.py")
runner = load_module("phase0_runner", "tests/regression/run_regression.py")


class Phase0FixtureTests(unittest.TestCase):
    def test_reference_and_read_generation_is_deterministic(self):
        sequence = fixture.deterministic_sequence("chr1", 50_000)
        self.assertEqual(len(sequence), 50_000)
        self.assertEqual(
            fixture.hashlib.sha256(sequence.encode()).hexdigest(),
            "4ad902181fb09cd40a8517963dac9b42046a2ae1dbf5293b6745257e175f2d6d",
        )

    def test_ligation_and_tag_positions_match_real_parser_contract(self):
        self.assertEqual(fixture.RNA_I1_TEMPLATE[15:23], "ACGTACGT")
        self.assertEqual(fixture.RNA_I1_TEMPLATE[53:61], "TGCATGCA")
        self.assertEqual(fixture.RNA_I1_TEMPLATE[91:99], "GATCGATC")
        self.assertEqual(fixture.DNA_SINGLE_I1_TEMPLATE[15:23], "ACGTACGT")
        self.assertEqual(fixture.DNA_SINGLE_I1_TEMPLATE[53:61], "TGCATGCA")
        self.assertEqual(fixture.DNA_SINGLE_I1_TEMPLATE[91:99], "GATCGATC")
        self.assertEqual(fixture.DNA_SINGLE_I2_TEMPLATE[14:18], "TACG")
        self.assertEqual(fixture.DNA_SINGLE_I2_TEMPLATE[18:26], "TACGTAAA")
        self.assertEqual(fixture.DNA_DUAL_I1_TEMPLATE[0:3], "AAA")
        self.assertEqual(fixture.DNA_DUAL_I1_TEMPLATE[3:11], "AGGCTATA")
        self.assertEqual(fixture.DNA_DUAL_I1_TEMPLATE[41:49], "ACGTACGT")
        self.assertEqual(fixture.DNA_DUAL_I1_TEMPLATE[79:87], "TGCATGCA")
        self.assertEqual(fixture.DNA_DUAL_I1_TEMPLATE[117:125], "GATCGATC")
        signature_lines = (
            REPO / "assets/dual_tag_artifact_23mers.fasta"
        ).read_text().splitlines()
        self.assertIn(fixture.DUAL_TAG_ARTIFACT, signature_lines)

    def test_baseline_provenance_is_exact(self):
        provenance = json.loads(
            (REPO / "tests/regression/contracts/v1.1.1/provenance.json").read_text()
        )
        self.assertEqual(provenance["release"], "v1.1.1")
        self.assertEqual(
            provenance["commit"], "40a383e80d945618952f8e5bfddb73ed3dc63af6"
        )
        self.assertEqual(
            provenance["tag_object"], "043fc873f9af89af49bf51203571d79a49ee1f01"
        )
        for scenario in ("rna_only", "dna_single", "dna_dual"):
            payload = json.loads(
                (REPO / f"tests/regression/contracts/v1.1.1/{scenario}.json").read_text()
            )
            self.assertEqual(payload["baseline_provenance"], provenance)

    def test_required_tool_versions_are_pinned(self):
        self.assertEqual(fixture.EXPECTED_STAR_VERSION, "2.7.11b")
        self.assertEqual(fixture.EXPECTED_BWA_MEM2_VERSION, "2.2.1")

    def test_runner_validates_all_golden_provenance(self):
        provenance = runner.validate_baseline_provenance(
            REPO / "tests/regression/contracts/v1.1.1"
        )
        self.assertEqual(provenance["commit"], runner.BASELINE_COMMIT)


class Phase0NormalizerTests(unittest.TestCase):
    def test_incomplete_contract_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete semantic regression contract"):
            normalizer.validate_contract(
                {
                    "scenario": "rna_only",
                    "directories": [],
                    "files": [],
                    "tables": {},
                    "gzip_text": {},
                    "bams": {},
                    "bigwigs": {},
                    "html": {},
                    "fastqc": {},
                    "multiqc": {},
                }
            )

    def test_paths_timestamps_and_gzip_headers_are_normalized(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            first_root = Path(first)
            second_root = Path(second)
            first_file = first_root / "records.tsv.gz"
            second_file = second_root / "records.tsv.gz"
            with gzip.GzipFile(first_file, "wb", mtime=1) as handle:
                handle.write(f"b\t{first_root}\t2026-01-01T01:02:03Z\na\n".encode())
            with gzip.GzipFile(second_file, "wb", mtime=2) as handle:
                handle.write(f"a\nb\t{second_root}\t2027-02-02T02:03:04Z\n".encode())
            self.assertNotEqual(first_file.read_bytes(), second_file.read_bytes())
            self.assertEqual(
                normalizer.semantic_gzip_text(first_file, (first_root,)),
                normalizer.semantic_gzip_text(second_file, (second_root,)),
            )

    def test_star_mapping_speed_inf_and_finite_values_are_runtime_metadata(self):
        relative_path = "rna_align/sample.Log.final.out"
        infinite = "       Mapping speed, Million of reads per hour |\tinf"
        finite = "       Mapping speed, Million of reads per hour |\t0.23"
        canonical_infinite = normalizer.normalize_scoped_runtime_line(
            relative_path, infinite
        )
        canonical_finite = normalizer.normalize_scoped_runtime_line(
            relative_path, finite
        )
        self.assertEqual(canonical_infinite, canonical_finite)
        self.assertIn(normalizer.STAR_MAPPING_SPEED_FIELD, canonical_infinite)
        self.assertTrue(canonical_infinite.endswith("<RUNTIME>"))

    def test_picard_started_on_weekday_and_timestamp_are_runtime_metadata(self):
        relative_path = "dna_align/sample.DuplicateMetrics.txt"
        legacy_expected = "# Started on: Wed <TIMESTAMP> CEST 2026"
        newly_observed = "# Started on: Thu Aug 27 10:27:22 CEST 2026"
        self.assertEqual(
            normalizer.normalize_scoped_runtime_line(relative_path, legacy_expected),
            normalizer.normalize_scoped_runtime_line(relative_path, newly_observed),
        )
        self.assertEqual(
            normalizer.normalize_scoped_runtime_line(relative_path, newly_observed),
            "# Started on: <TIMESTAMP>",
        )

    def test_expected_and_observed_contracts_are_normalized_symmetrically(self):
        expected = {
            "text": {
                "rna_align/sample.Log.final.out": [
                    "       Mapping speed, Million of reads per hour |\tinf"
                ],
                "dna_align/sample.DuplicateMetrics.txt": [
                    "# Started on: Wed <TIMESTAMP> CEST 2026"
                ],
            }
        }
        observed = {
            "text": {
                "rna_align/sample.Log.final.out": [
                    "       Mapping speed, Million of reads per hour |\t0.23"
                ],
                "dna_align/sample.DuplicateMetrics.txt": [
                    "# Started on: Thu Aug 27 10:27:22 CEST 2026"
                ],
            }
        }
        self.assertEqual(
            normalizer.canonicalize_contract_runtime_metadata(expected),
            normalizer.canonicalize_contract_runtime_metadata(observed),
        )

    def test_retired_process_host_dependencies_are_normalized_symmetrically(self):
        expected = {
            "tables": {
                normalizer.RUNTIME_CONTRACT_PATH: {
                    "delimiter": "tab",
                    "rows": [
                        ["tool", "configured_path", "exists", "currently_used"],
                        ["python3", "<PATH>", "true", "yes"],
                        ["cutadapt", "<PATH>", "true", "yes"],
                        ["trim_galore", "<PATH>", "true", "yes"],
                        ["codon", "<PATH>", "true", "yes"],
                        ["pigz", "<PATH>", "true", "yes"],
                        ["STAR", "<PATH>", "true", "yes"],
                        ["bedGraphToBigWig", "<PATH>", "true", "yes"],
                        ["[runtime_environment]"],
                        ["runtime_env_prefix", "<PATH>"],
                        ["codon_home", "<PATH>"],
                        ["runtime_tmpdir", "<PATH>"],
                        ["[host_codon_seq_preflight]"],
                        ["required_codon_version=0.16.3"],
                        ["required_seq_version=0.11.3"],
                    ],
                }
            }
        }
        observed = {
            "tables": {
                normalizer.RUNTIME_CONTRACT_PATH: {
                    "delimiter": "tab",
                    "rows": [
                        ["tool", "configured_path", "exists", "currently_used"],
                        ["python3", "<PATH>", "true", "yes"],
                        ["[runtime_environment]"],
                        ["runtime_env_prefix", "<PATH>"],
                        ["runtime_tmpdir", "<PATH>"],
                    ],
                }
            }
        }
        self.assertEqual(
            normalizer.canonicalize_contract_runtime_metadata(expected),
            normalizer.canonicalize_contract_runtime_metadata(observed),
        )

    def test_nonretired_runtime_contract_changes_remain_visible(self):
        expected = {
            "tables": {
                normalizer.RUNTIME_CONTRACT_PATH: {
                    "delimiter": "tab",
                    "rows": [["samtools", "<PATH>", "true", "yes"]],
                }
            }
        }
        observed = {
            "tables": {
                normalizer.RUNTIME_CONTRACT_PATH: {
                    "delimiter": "tab",
                    "rows": [["samtools", "<PATH>", "false", "yes"]],
                }
            }
        }
        self.assertNotEqual(
            normalizer.canonicalize_contract_runtime_metadata(expected),
            normalizer.canonicalize_contract_runtime_metadata(observed),
        )

    def test_star_and_picard_biological_metric_differences_remain_visible(self):
        stable_runtime_lines = {
            "rna_align/sample.Log.final.out": [
                "       Mapping speed, Million of reads per hour |\tinf"
            ],
            "dna_align/sample.DuplicateMetrics.txt": [
                "# Started on: Wed <TIMESTAMP> CEST 2026"
            ],
        }
        expected = {
            "text": stable_runtime_lines
            | {
                "rna_align/sample.Log.final.out": [
                    *stable_runtime_lines["rna_align/sample.Log.final.out"],
                    "                          Number of input reads |\t64",
                ],
                "dna_align/sample.DuplicateMetrics.txt": [
                    *stable_runtime_lines["dna_align/sample.DuplicateMetrics.txt"],
                    "PHASE0_SYNTHETIC\t0\t64\t0\t0\t0\t63\t63\t0.984375",
                ],
            }
        }
        star_difference = {
            "text": expected["text"]
            | {
                "rna_align/sample.Log.final.out": [
                    *stable_runtime_lines["rna_align/sample.Log.final.out"],
                    "                          Number of input reads |\t65",
                ]
            }
        }
        picard_difference = {
            "text": expected["text"]
            | {
                "dna_align/sample.DuplicateMetrics.txt": [
                    *stable_runtime_lines["dna_align/sample.DuplicateMetrics.txt"],
                    "PHASE0_SYNTHETIC\t0\t64\t0\t0\t0\t62\t62\t0.968750",
                ]
            }
        }
        canonical_expected = normalizer.canonicalize_contract_runtime_metadata(expected)
        self.assertNotEqual(
            canonical_expected,
            normalizer.canonicalize_contract_runtime_metadata(star_difference),
        )
        self.assertNotEqual(
            canonical_expected,
            normalizer.canonicalize_contract_runtime_metadata(picard_difference),
        )

    def test_html_comparison_uses_visible_data_not_script_or_style(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.html"
            path.write_text(
                "<html><style>.x{color:red}</style><body><h1>Report</h1>"
                "<p>Generated 2026-08-26T12:00:00+00:00</p>"
                "<table><tr><td>reads</td><td>64</td></tr></table>"
                "<script>volatile='ignore me'</script></body></html>"
            )
            self.assertEqual(
                normalizer.semantic_html(path, (Path(directory),)),
                ["Report", "Generated <TIMESTAMP>", "reads", "64"],
            )

    @unittest.skipUnless(shutil.which("samtools"), "samtools is required")
    def test_bam_records_headers_flags_coordinates_cigar_and_tags_are_semantic(self):
        samtools = Path(shutil.which("samtools"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sam = root / "input.sam"
            bam = root / "input.bam"
            sam.write_text(
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@PG\tID:run123\tPN:samtools\tVN:1.23\tCL:/tmp/run123\n"
                "read1\t99\tchr1\t11\t60\t4M\t=\t21\t14\tACGT\tIIII\tCB:Z:CELL\tNM:i:0\n"
            )
            subprocess.run(
                [str(samtools), "view", "-b", "-o", str(bam), str(sam)], check=True
            )
            decoded = normalizer.semantic_bam(bam, samtools, (root,))
            self.assertIn("@HD\tSO:coordinate", decoded["header"])
            self.assertIn("@PG\tPN:samtools", decoded["header"])
            self.assertEqual(decoded["records"][0][0:8], ["read1", "99", "chr1", "11", "60", "4M", "=", "21"])
            self.assertEqual(decoded["records"][0][-2:], ["CB:Z:CELL", "NM:i:0"])

    def test_star_bam_header_thread_allocation_is_runtime_metadata(self):
        first = normalizer.canonical_header_line(
            "@CO\tuser command line: STAR --runThreadN 2 --soloFeatures GeneFull",
            (),
        )
        second = normalizer.canonical_header_line(
            "@CO\tuser command line: STAR --runThreadN 16 --soloFeatures GeneFull",
            (),
        )
        self.assertEqual(first, second)
        self.assertIn("--runThreadN <THREADS>", first)

    def test_star_bam_header_staged_and_absolute_genome_dirs_are_runtime_paths(self):
        staged = normalizer.canonical_header_line(
            "@CO\tuser command line: STAR --genomeDir star --soloFeatures GeneFull",
            (),
        )
        absolute = normalizer.canonical_header_line(
            "@CO\tuser command line: STAR --genomeDir /refs/hg38/star --soloFeatures GeneFull",
            (),
        )
        self.assertEqual(staged, absolute)
        self.assertIn("--genomeDir <PATH>", staged)

    def test_bigwig_intervals_and_values_are_semantic(self):
        try:
            import pyBigWig
        except ImportError:
            self.skipTest("pyBigWig is required")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coverage.bw"
            with pyBigWig.open(str(path), "w") as bigwig:
                bigwig.addHeader([("chr1", 100)])
                bigwig.addEntries(["chr1", "chr1"], [10, 20], ends=[15, 30], values=[1.5, 2.25])
            self.assertEqual(
                normalizer.semantic_bigwig(path),
                {
                    "chromosomes": {"chr1": 100},
                    "intervals": {"chr1": [[10, 15, 1.5], [20, 30, 2.25]]},
                },
            )


if __name__ == "__main__":
    unittest.main()
