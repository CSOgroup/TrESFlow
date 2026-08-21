import importlib.util
import gzip
import importlib
import os
import tempfile
import unittest
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "bin"))


def load_module(name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUN_TAG = load_module("run_tag", "bin/run_tag.py")
RUN_TAG_UMI = load_module("run_tag_umi", "bin/run_tag_umi.py")
RUN_TAG_LIG3 = load_module("run_tag_lig3", "bin/run_tag_lig3.py")
SPLIT_RNA = load_module("run_split_reads_rna", "bin/run_split_reads_rna.py")
SPLIT_DNA = load_module("run_split_reads_dna", "bin/run_split_reads_dna.py")
FQ_TO_SAM = load_module("run_fq_to_sam", "bin/run_fq_to_sam.py")
RUN_TRIM_GALORE = load_module("run_trim_galore", "bin/run_trim_galore.py")
UTILS = importlib.import_module("tresflow_fastq_utils")


class FastqCompressionTests(unittest.TestCase):
    def test_mock_trim_galore_outputs_uncompressed_fastqs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            r1 = tmp_path / "source_R1.fastq.gz"
            r2 = tmp_path / "source_R2.fastq.gz"
            output_r1 = tmp_path / "trimmed_R1.fq"
            output_r2 = tmp_path / "trimmed_R2.fq"
            record = b"@r1\nACGT\n+\n!!!!\n"

            for source in (r1, r2):
                with gzip.open(source, "wb") as handle:
                    handle.write(record)

            RUN_TRIM_GALORE.mock_trim(
                SimpleNamespace(r1=r1, r2=r2, output_r1=output_r1, output_r2=output_r2)
            )

            self.assertEqual(output_r1.read_bytes(), record)
            self.assertEqual(output_r2.read_bytes(), record)
            self.assertNotEqual(output_r1.read_bytes()[:2], b"\x1f\x8b")
            self.assertNotEqual(output_r2.read_bytes()[:2], b"\x1f\x8b")

    def test_real_trim_galore_explicitly_disables_gzip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            trim_galore = tmp_path / "trim_galore"
            trim_galore.write_text("#!/bin/sh\n", encoding="utf-8")
            trim_galore.chmod(0o755)
            r1 = tmp_path / "source_R1.fastq.gz"
            r2 = tmp_path / "source_R2.fastq.gz"
            output_r1 = tmp_path / "trimmed_R1.fq"
            output_r2 = tmp_path / "trimmed_R2.fq"
            calls = []

            def fake_run(command, check):
                calls.append((command, check))
                output_dir = Path(command[command.index("--output_dir") + 1])
                (output_dir / "source_R1_val_1.fq").write_bytes(b"r1")
                (output_dir / "source_R2_val_2.fq").write_bytes(b"r2")

            args = SimpleNamespace(
                r1=r1,
                r2=r2,
                output_r1=output_r1,
                output_r2=output_r2,
                trim_galore_bin=trim_galore,
                quality=10,
                cores=2,
                length=20,
            )

            with mock.patch.dict(os.environ, {"TMPDIR": str(tmp_path)}), mock.patch.object(
                RUN_TRIM_GALORE.subprocess, "run", side_effect=fake_run
            ):
                RUN_TRIM_GALORE.real_trim(args)

            command, check = calls[0]
            self.assertTrue(check)
            self.assertIn("--dont_gzip", command)
            self.assertNotIn("--gzip", command)
            self.assertEqual(output_r1.read_bytes(), b"r1")
            self.assertEqual(output_r2.read_bytes(), b"r2")

    def test_tag_wrappers_fail_on_compression_mismatch(self):
        for module in (RUN_TAG, RUN_TAG_UMI, RUN_TAG_LIG3):
            with self.subTest(module=module.__name__):
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmp_path = Path(tmpdir)
                    source = tmp_path / "source.fastq"
                    destination = tmp_path / "destination.fastq.gz"
                    source.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")

                    with self.assertRaisesRegex(
                        RuntimeError,
                        "Python wrappers must not recompress production FASTQs",
                    ) as context:
                        module.strict_move_fastq(source, destination)

                    message = str(context.exception)
                    self.assertIn(str(source), message)
                    self.assertIn(str(destination), message)
                    self.assertTrue(source.exists())
                    self.assertFalse(destination.exists())

    def test_split_fastq_names_normalize_without_changing_compression_state(self):
        cases = {
            "sample_Normal_R1.fq": "sample_Normal_R1.fastq",
            "sample_Normal_R2.fq.gz": "sample_Normal_R2.fastq.gz",
            "sample_Normal_R1.fastq": "sample_Normal_R1.fastq",
            "sample_Normal_R2.fastq.gz": "sample_Normal_R2.fastq.gz",
        }
        for source, expected in cases.items():
            self.assertEqual(UTILS.normalize_split_fastq_name(source), expected)

    def test_split_publication_compression_is_a_separate_non_destructive_process(self):
        module_text = (REPO_ROOT / "modules/local/compress_split_fastqs/main.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn('-c -p "${task.cpus}"', module_text)
        self.assertIn('> "\\$(basename "\\$fastq").gz"', module_text)
        self.assertNotIn('pigz -f', module_text)
        self.assertIn('.fastq.gz', module_text)

    def test_canonical_cell_id_drops_modality_specific_sb_without_replacing_technical_cb(self):
        cell_barcode = "ACGTACGTTGCATGCAGATCGATC"
        rna_comment = f"CB:Z:CAGT{cell_barcode}\tRG:Z:CAGT{cell_barcode}\tUM:Z:TTTT\tSB:Z:CAGT"
        dna_comment = f"CB:Z:AAA{cell_barcode}\tRG:Z:AAA{cell_barcode}\tMO:Z:AGGCTATA\tSB:Z:AAA"

        rna_canonical = SPLIT_RNA.canonicalize_fastq_comment("sample1", "Normal", rna_comment)
        dna_canonical = SPLIT_DNA.canonicalize_dna_fastq_comment(
            "sample1",
            "Normal",
            "AV240401:AVT0507:2528453125:1:11104:5031:3419:ACGT",
            dna_comment,
        )
        expected = f"sample1_Normal_{cell_barcode}"

        self.assertIn(f"CB:Z:{cell_barcode}", rna_canonical)
        self.assertIn(f"RG:Z:{cell_barcode}", rna_canonical)
        self.assertIn(f"XI:Z:{expected}", rna_canonical)
        self.assertIn(f"CB:Z:{cell_barcode}", dna_canonical)
        self.assertIn("RG:Z:AV240401:AVT0507:2528453125:L1", dna_canonical)
        self.assertIn(f"XI:Z:{expected}", dna_canonical)
        self.assertIn("SB:Z:CAGT", rna_canonical)
        self.assertIn("SB:Z:AAA", dna_canonical)

    def test_fq_to_sam_uses_raw_technical_cb_for_star_cb_length(self):
        cell_barcode = "ACGTACGTTGCATGCAGATCGATC"
        umi = "TTTTGGGGAA"
        canonical = f"Isa_VeryLongGroupName_{cell_barcode}"
        comment = f"CB:Z:{cell_barcode}\tRG:Z:{cell_barcode}\tUM:Z:{umi}\tSB:Z:CAGT\tXI:Z:{canonical}"

        cr_value, other_tags = FQ_TO_SAM.extract_cr_and_others(comment)

        self.assertEqual(cr_value, cell_barcode + umi)
        self.assertEqual(len(cr_value) - len(umi), len(cell_barcode))
        self.assertNotIn("Isa_VeryLongGroupName", cr_value)
        self.assertIn(f"XI:Z:{canonical}", other_tags)


if __name__ == "__main__":
    unittest.main()
