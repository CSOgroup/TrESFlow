import importlib.util
import importlib
import tempfile
import unittest
import sys
from pathlib import Path


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
UTILS = importlib.import_module("tresflow_fastq_utils")


class FastqCompressionTests(unittest.TestCase):
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

    def test_split_fastq_names_normalize_to_fastq_gz(self):
        cases = {
            "sample_Normal_R1.fq": "sample_Normal_R1.fastq",
            "sample_Normal_R2.fq.gz": "sample_Normal_R2.fastq.gz",
            "sample_Normal_R1.fastq": "sample_Normal_R1.fastq",
            "sample_Normal_R2.fastq.gz": "sample_Normal_R2.fastq.gz",
        }
        for source, expected in cases.items():
            self.assertEqual(SPLIT_RNA.normalize_split_fastq_name(source), expected)
            self.assertEqual(SPLIT_DNA.normalize_split_fastq_name(source), expected)

    def test_final_compression_uses_pigz_thread_count(self):
        for module in (SPLIT_RNA, SPLIT_DNA):
            with self.subTest(module=module.__name__):
                with tempfile.TemporaryDirectory() as tmpdir:
                    source = Path(tmpdir) / "sample_Normal_R1.fastq"
                    source.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
                    calls = []
                    original_run = UTILS.subprocess.run

                    def fake_run(command, check):
                        calls.append((command, check))
                        Path(str(source) + ".gz").write_bytes(b"compressed")
                        source.unlink()

                    try:
                        UTILS.subprocess.run = fake_run
                        module.compress_fastq_with_pigz(source, 6, "/usr/bin/pigz")
                    finally:
                        UTILS.subprocess.run = original_run

                    self.assertEqual(calls, [(["/usr/bin/pigz", "-p", "6", str(source)], True)])

    def test_canonical_cell_id_drops_modality_specific_sb_without_replacing_technical_cb(self):
        cell_barcode = "ACGTACGTTGCATGCAGATCGATC"
        rna_comment = f"CB:Z:CAGT{cell_barcode}\tRG:Z:CAGT{cell_barcode}\tUM:Z:TTTT\tSB:Z:CAGT"
        dna_comment = f"CB:Z:AAA{cell_barcode}\tRG:Z:AAA{cell_barcode}\tMO:Z:AGGCTATA\tSB:Z:AAA"

        rna_canonical = SPLIT_RNA.canonicalize_fastq_comment("sample1", "Normal", rna_comment)
        dna_canonical = SPLIT_DNA.canonicalize_fastq_comment("sample1", "Normal", dna_comment)
        expected = f"sample1_Normal_{cell_barcode}"

        self.assertIn(f"CB:Z:{cell_barcode}", rna_canonical)
        self.assertIn(f"RG:Z:{cell_barcode}", rna_canonical)
        self.assertIn(f"XI:Z:{expected}", rna_canonical)
        self.assertIn(f"CB:Z:{cell_barcode}", dna_canonical)
        self.assertIn(f"RG:Z:{cell_barcode}", dna_canonical)
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
