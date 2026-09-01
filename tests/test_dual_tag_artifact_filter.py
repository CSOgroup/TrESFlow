import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
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


def read_fastq_records(path):
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    if len(lines) % 4:
        raise AssertionError(f"Malformed FASTQ fixture output: {path}")
    return ["".join(lines[index : index + 4]) for index in range(0, len(lines), 4)]


def fastq_record_name(record):
    return (
        record.splitlines()[0]
        .split()[0]
        .removeprefix("@")
        .removesuffix("/1")
        .removesuffix("/2")
    )


def fastq_record_sequence(record):
    return record.splitlines()[1]


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
        rna_workflow = (REPO_ROOT / "subworkflows/local/rna_core/main.nf").read_text(
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
        self.assertIn("TRIM_DNA_FASTQS.out.trimmed.branch", dna_workflow)
        self.assertNotIn("TAG_DNA_CELL_BARCODE.out.tagged.branch", dna_workflow)
        self.assertIn(".mix(DUAL_TAG_ARTIFACT_FILTER.out.filtered)", dna_workflow)
        self.assertIn(".join(ch_trimmed_for_split)", dna_workflow)
        self.assertNotIn(".join(TRIM_DNA_FASTQS.out.trimmed)", dna_workflow)
        self.assertLess(
            dna_workflow.index(
                "TRIM_DNA_FASTQS(TAG_DNA_CELL_BARCODE.out.tagged, trimHelperScript)"
            ),
            dna_workflow.index("TRIM_DNA_FASTQS.out.trimmed.branch"),
        )
        self.assertLess(
            dna_workflow.index("TRIM_DNA_FASTQS.out.trimmed.branch"),
            dna_workflow.index("DUAL_TAG_ARTIFACT_FILTER("),
        )
        self.assertNotIn("DUAL_TAG_ARTIFACT_FILTER", rna_workflow)

    def test_cutadapt_is_a_portable_process_dependency(self):
        runtime_module = (REPO_ROOT / "modules/local/runtime_support/main.nf").read_text(
            encoding="utf-8"
        )
        runtime_support = (REPO_ROOT / "lib/RuntimeSupport.groovy").read_text(
            encoding="utf-8"
        )
        process_module = (
            REPO_ROOT / "modules/local/dual_tag_artifact_filter/main.nf"
        ).read_text(encoding="utf-8")

        self.assertNotIn("CUTADAPT_BIN", runtime_module)
        self.assertNotIn("[name: 'cutadapt', binary: 'cutadapt']", runtime_support)
        self.assertIn("environment-cutadapt.yml", process_module)
        self.assertIn("quay.io/biocontainers/cutadapt@sha256:", process_module)
        self.assertIn('path helperScript, stageAs:', process_module)
        self.assertIn('python3 "tresflow/bin/run_dual_tag_artifact_filter.py"', process_module)
        self.assertNotIn('--cutadapt-bin', process_module)

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
        cls.trim_galore = shutil.which("trim_galore")

    def test_real_trim_galore_then_filter_salvages_prefix_and_rejects_only_residuals(self):
        if not self.trim_galore or not self.cutadapt:
            self.skipTest(
                "trim_galore and cutadapt are required for the direct wrapper-chain test"
            )

        signatures = FILTER.load_and_validate_signatures(ASSET)
        signature = signatures[0]
        one_mismatch = (
            signature[:11]
            + ("T" if signature[11] != "T" else "A")
            + signature[12:]
        )
        genomic_prefix = "GCTACCGTTCAGGATCCGTACTGGCATGCC"
        illumina_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        neutral = "CGTTCGAGTCCGATGCTACCGTTCGAGTCCGATGCTACCGTTC"

        r1_sequences = [
            genomic_prefix + illumina_adapter + signature,
            "GCGC" + signature + "TCTC",
            neutral,
            signature[:22],
            neutral,
            neutral,
        ]
        r2_sequences = [
            neutral,
            neutral,
            "GCGC" + signature + "TCTC",
            neutral,
            "GCGC" + one_mismatch + "TCTC",
            neutral,
        ]
        names = [
            "salvage",
            "residual_r1",
            "residual_r2",
            "partial22",
            "mismatch",
            "clean",
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_r1 = root / "input-R1.fastq"
            input_r2 = root / "input-R2.fastq"
            trimmed_r1 = root / "trimmed-R1.fq"
            trimmed_r2 = root / "trimmed-R2.fq"
            filtered_r1 = root / "filtered-R1.fastq"
            filtered_r2 = root / "filtered-R2.fastq"
            output_json = root / "sample.dual_tag_artifact_filter.cutadapt.json"
            output_summary = root / "sample.dual_tag_artifact_filter.summary.tsv"
            input_r1.write_text(
                "".join(
                    fastq_record(f"{name}/1", sequence)
                    for name, sequence in zip(names, r1_sequences)
                ),
                encoding="utf-8",
            )
            input_r2.write_text(
                "".join(
                    fastq_record(f"{name}/2", sequence)
                    for name, sequence in zip(names, r2_sequences)
                ),
                encoding="utf-8",
            )

            trim_tmpdir = root / "trim_tmp"
            trim_tmpdir.mkdir()
            trim_environment = os.environ.copy()
            trim_environment["TMPDIR"] = str(trim_tmpdir)
            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "bin/run_trim_galore.py"),
                    "--mode",
                    "real",
                    "--r1",
                    str(input_r1),
                    "--r2",
                    str(input_r2),
                    "--trim-galore-bin",
                    self.trim_galore,
                    "--quality",
                    "10",
                    "--cores",
                    "1",
                    "--length",
                    "20",
                    "--output-r1",
                    str(trimmed_r1),
                    "--output-r2",
                    str(trimmed_r2),
                ],
                check=True,
                capture_output=True,
                text=True,
                env=trim_environment,
            )

            trimmed_r1_records = read_fastq_records(trimmed_r1)
            trimmed_r2_records = read_fastq_records(trimmed_r2)
            self.assertEqual(len(trimmed_r1_records), 6)
            self.assertEqual(len(trimmed_r2_records), 6)
            self.assertEqual(fastq_record_sequence(trimmed_r1_records[0]), genomic_prefix)
            self.assertNotIn(signature, fastq_record_sequence(trimmed_r1_records[0]))

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
                    str(trimmed_r1),
                    "--r2",
                    str(trimmed_r2),
                    "--signature-fasta",
                    str(ASSET),
                    "--cutadapt-bin",
                    self.cutadapt,
                    "--cpus",
                    "1",
                    "--output-r1",
                    str(filtered_r1),
                    "--output-r2",
                    str(filtered_r2),
                    "--output-json",
                    str(output_json),
                    "--output-summary",
                    str(output_summary),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            rejected_names = {"residual_r1", "residual_r2"}
            expected_r1 = "".join(
                record
                for record in trimmed_r1_records
                if fastq_record_name(record) not in rejected_names
            )
            expected_r2 = "".join(
                record
                for record in trimmed_r2_records
                if fastq_record_name(record) not in rejected_names
            )
            filtered_r1_records = read_fastq_records(filtered_r1)
            filtered_r2_records = read_fastq_records(filtered_r2)
            summary_lines = output_summary.read_text(encoding="utf-8").splitlines()
            summary = dict(zip(summary_lines[0].split("\t"), summary_lines[1].split("\t")))

            self.assertEqual(filtered_r1.read_text(encoding="utf-8"), expected_r1)
            self.assertEqual(filtered_r2.read_text(encoding="utf-8"), expected_r2)
            self.assertEqual(
                [fastq_record_name(record) for record in filtered_r1_records],
                ["salvage", "partial22", "mismatch", "clean"],
            )
            self.assertEqual(
                [fastq_record_name(record) for record in filtered_r1_records],
                [fastq_record_name(record) for record in filtered_r2_records],
            )
            self.assertEqual(int(summary["input_pairs"]), 6)
            self.assertEqual(int(summary["retained_pairs"]), 4)
            self.assertEqual(int(summary["rejected_pairs"]), 2)
            self.assertEqual(
                int(summary["input_pairs"]),
                int(summary["retained_pairs"]) + int(summary["rejected_pairs"]),
            )

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
