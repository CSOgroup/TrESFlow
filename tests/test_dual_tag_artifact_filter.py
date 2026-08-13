import argparse
import importlib.util
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
ASSET = REPO_ROOT / "assets/dual_tag_artifact_23mers.fasta"
EXPECTED_SHA256 = "67c6f1789ef5e36492562203ac38fc13fa901058047ed2bd37b304d85a30ae0f"


def load_wrapper():
    spec = importlib.util.spec_from_file_location(
        "dual_tag_artifact_filter", REPO_ROOT / "bin/run_dual_tag_artifact_filter.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FILTER = load_wrapper()


def fastq_record(name, sequence):
    return f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"


class DualTagArtifactFilterUnitTests(unittest.TestCase):
    def test_runtime_asset_is_the_audited_signature_set(self):
        signatures = FILTER.load_and_validate_signatures(ASSET)
        self.assertEqual(FILTER.file_sha256(ASSET), EXPECTED_SHA256)
        self.assertEqual(len(signatures), 48)
        self.assertEqual(len(set(signatures)), 48)
        self.assertTrue(all(len(sequence) == 23 for sequence in signatures))
        self.assertTrue(all(set(sequence) <= set("ACGT") for sequence in signatures))

    def test_cutadapt_command_has_exact_detection_and_discard_semantics(self):
        args = argparse.Namespace(
            cutadapt_bin=Path("/bin/true"),
            cpus=4,
            signature_fasta=Path("signatures.fasta"),
            output_json=Path("sample.json"),
            output_r1=Path("clean-R1.fastq"),
            output_r2=Path("clean-R2.fastq"),
            r1=Path("input-R1.fastq"),
            r2=Path("input-R2.fastq"),
        )
        self.assertEqual(
            FILTER.build_cutadapt_command(args),
            [
                "/bin/true",
                "-j",
                "4",
                "-e",
                "0",
                "-O",
                "23",
                "--no-indels",
                "--action=none",
                "-b",
                "file:signatures.fasta",
                "-B",
                "file:signatures.fasta",
                "--discard-trimmed",
                "--pair-filter=any",
                "--json",
                "sample.json",
                "-o",
                "clean-R1.fastq",
                "-p",
                "clean-R2.fastq",
                "input-R1.fastq",
                "input-R2.fastq",
            ],
        )

    def test_parameter_defaults_and_dual_only_routing(self):
        config = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")
        schema = json.loads((REPO_ROOT / "nextflow_schema.json").read_text(encoding="utf-8"))
        dna_workflow = (REPO_ROOT / "subworkflows/local/dna_core/main.nf").read_text(
            encoding="utf-8"
        )
        parameter = schema["$defs"]["execution_options"]["properties"][
            "filter_dual_tag_artifacts"
        ]

        self.assertIn("filter_dual_tag_artifacts = true", config)
        self.assertEqual(parameter, {
            "type": "boolean",
            "description": parameter["description"],
            "default": True,
        })
        self.assertIn(
            "filter: params.filter_dual_tag_artifacts && meta.dna_tagmentation == 'dual'",
            dna_workflow,
        )
        self.assertIn("bypass: true", dna_workflow)
        self.assertEqual(dna_workflow.count("DUAL_TAG_ARTIFACT_FILTER("), 1)
        self.assertIn(".mix(DUAL_TAG_ARTIFACT_FILTER.out.filtered)", dna_workflow)

    def test_cutadapt_is_an_explicit_runtime_dependency(self):
        runtime_module = (REPO_ROOT / "modules/local/runtime_support/main.nf").read_text(
            encoding="utf-8"
        )
        runtime_support = (REPO_ROOT / "lib/RuntimeSupport.groovy").read_text(
            encoding="utf-8"
        )
        process_module = (
            REPO_ROOT / "modules/local/dual_tag_artifact_filter/main.nf"
        ).read_text(encoding="utf-8")

        self.assertIn("CUTADAPT_BIN", runtime_module)
        self.assertIn("[name: 'cutadapt', binary: 'cutadapt']", runtime_support)
        self.assertIn('--cutadapt-bin "\\$CUTADAPT_BIN"', process_module)

    def test_zero_input_mock_has_safe_fraction_and_balanced_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_r1 = root / "input-R1.fastq"
            input_r2 = root / "input-R2.fastq"
            input_r1.write_text("", encoding="utf-8")
            input_r2.write_text("", encoding="utf-8")
            args = argparse.Namespace(
                mode="mock",
                sample_id="empty",
                tagmentation="dual",
                r1=input_r1,
                r2=input_r2,
                signature_fasta=ASSET,
                cutadapt_bin=None,
                cpus=1,
                output_r1=root / "clean-R1.fastq",
                output_r2=root / "clean-R2.fastq",
                output_json=root / "empty.json",
                output_summary=root / "empty.tsv",
            )
            FILTER.run(args)
            lines = args.output_summary.read_text(encoding="utf-8").splitlines()
            summary = dict(zip(lines[0].split("\t"), lines[1].split("\t")))

        self.assertEqual(summary["input_pairs"], "0")
        self.assertEqual(summary["retained_pairs"], "0")
        self.assertEqual(summary["rejected_pairs"], "0")
        self.assertEqual(summary["rejected_fraction"], "0.00000000")

    def test_alignment_command_is_unchanged(self):
        alignment = (REPO_ROOT / "scripts/core_runtime/AlignDNA.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            '"${BWA_MEM2_BIN}" mem -t ${threads} -C -o ${RGID}_TEMP.sam '
            '${path_bwarefDB} ${R1} ${R2}',
            alignment,
        )
        self.assertNotIn(" -k ", alignment)


class DualTagArtifactFilterCutadaptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cutadapt = shutil.which("cutadapt")

    def test_exact_signatures_discard_pairs_and_retained_records_are_unchanged(self):
        if not self.cutadapt:
            self.skipTest("cutadapt is not available in PATH")

        signatures = FILTER.load_and_validate_signatures(ASSET)
        signature = signatures[0]
        one_mismatch = signature[:11] + ("T" if signature[11] != "T" else "A") + signature[12:]
        self.assertNotIn(one_mismatch, signatures)
        neutral = "GATTACAGATTACAGATTACAGATTACAGAT"

        r1_records = [
            fastq_record("r1_only_5prime/1", signature + "AAA"),
            fastq_record("r2_only_internal/1", neutral),
            fastq_record("both_r1_3prime/1", "AAA" + signature),
            fastq_record("only_22/1", signature[:22]),
            fastq_record("one_mismatch/1", one_mismatch),
            fastq_record("no_signature/1", neutral),
        ]
        r2_records = [
            fastq_record("r1_only_5prime/2", neutral),
            fastq_record("r2_only_internal/2", "AAA" + signature + "CCC"),
            fastq_record("both_r1_3prime/2", signature + "CCC"),
            fastq_record("only_22/2", neutral),
            fastq_record("one_mismatch/2", neutral),
            fastq_record("no_signature/2", neutral),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_r1 = root / "input-R1.fastq"
            input_r2 = root / "input-R2.fastq"
            output_r1 = root / "clean-R1.fastq"
            output_r2 = root / "clean-R2.fastq"
            output_json = root / "sample.dual_tag_artifact_filter.cutadapt.json"
            output_summary = root / "sample.dual_tag_artifact_filter.summary.tsv"
            input_r1.write_text("".join(r1_records), encoding="utf-8")
            input_r2.write_text("".join(r2_records), encoding="utf-8")

            subprocess.run(
                [
                    str(REPO_ROOT / "bin/run_dual_tag_artifact_filter.py"),
                    "--mode",
                    "real",
                    "--sample-id",
                    "sample",
                    "--tagmentation",
                    "dual",
                    "--r1",
                    str(input_r1),
                    "--r2",
                    str(input_r2),
                    "--signature-fasta",
                    str(ASSET),
                    "--cutadapt-bin",
                    self.cutadapt,
                    "--cpus",
                    "1",
                    "--output-r1",
                    str(output_r1),
                    "--output-r2",
                    str(output_r2),
                    "--output-json",
                    str(output_json),
                    "--output-summary",
                    str(output_summary),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            expected_r1 = "".join(r1_records[3:])
            expected_r2 = "".join(r2_records[3:])
            report = json.loads(output_json.read_text(encoding="utf-8"))
            summary_lines = output_summary.read_text(encoding="utf-8").splitlines()
            summary = dict(zip(summary_lines[0].split("\t"), summary_lines[1].split("\t")))

            self.assertEqual(output_r1.read_text(encoding="utf-8"), expected_r1)
            self.assertEqual(output_r2.read_text(encoding="utf-8"), expected_r2)
            self.assertEqual(report["read_counts"]["input"], 6)
            self.assertEqual(report["read_counts"]["output"], 3)
            self.assertEqual(report["read_counts"]["filtered"]["discard_trimmed"], 3)
            self.assertEqual(report["read_counts"]["read1_with_adapter"], 2)
            self.assertEqual(report["read_counts"]["read2_with_adapter"], 2)
            self.assertEqual(int(summary["input_pairs"]), 6)
            self.assertEqual(int(summary["retained_pairs"]), 3)
            self.assertEqual(int(summary["rejected_pairs"]), 3)
            self.assertEqual(int(summary["input_pairs"]), int(summary["retained_pairs"]) + int(summary["rejected_pairs"]))
            self.assertEqual(int(summary["r1_with_signature"]), 2)
            self.assertEqual(int(summary["r2_with_signature"]), 2)
            self.assertEqual(summary["signature_fasta_sha256"], EXPECTED_SHA256)


if __name__ == "__main__":
    unittest.main()
