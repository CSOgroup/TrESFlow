import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
QC_SCRIPT = REPO_ROOT / "scripts/core_runtime/SamtoolsBamQc.sh"


class SamtoolsBamQcTests(unittest.TestCase):
    def setUp(self):
        self.samtools = shutil.which("samtools")
        if not self.samtools:
            self.skipTest("samtools is not installed")

    def make_bam(self, root):
        source_sam = root / "source.sam"
        source_sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "pair1\t99\tchr1\t101\t60\t20M\t=\t151\t70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
            "pair1\t147\tchr1\t151\t60\t20M\t=\t101\t-70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n",
            encoding="utf-8",
        )
        bam = root / "synthetic.bam"
        subprocess.run(
            [self.samtools, "view", "--bam", "--output", str(bam), str(source_sam)],
            check=True,
        )
        subprocess.run([self.samtools, "index", str(bam)], check=True)
        return bam

    def run_qc(self, bam, prefix, run_idxstats):
        env = os.environ.copy()
        env["SAMTOOLS_BIN"] = self.samtools
        return subprocess.run(
            [
                "bash",
                str(QC_SCRIPT),
                str(bam),
                prefix.name,
                str(run_idxstats).lower(),
                "2",
            ],
            cwd=prefix.parent,
            env=env,
            capture_output=True,
            text=True,
        )

    def test_combined_qc_matches_individual_samtools_commands(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam = self.make_bam(root)
            prefix = root / "dna.synthetic"

            result = self.run_qc(bam, prefix, True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (root / "dna.synthetic.flagstat").read_text(encoding="utf-8"),
                subprocess.run(
                    [self.samtools, "flagstat", "--threads", "2", str(bam)],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout,
            )
            self.assertEqual(
                (root / "dna.synthetic.stats").read_text(encoding="utf-8"),
                subprocess.run(
                    [self.samtools, "stats", "--threads", "2", str(bam)],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout,
            )
            self.assertEqual(
                (root / "dna.synthetic.idxstats").read_text(encoding="utf-8"),
                subprocess.run(
                    [self.samtools, "idxstats", "--threads", "1", str(bam)],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout,
            )
            self.assertEqual(
                (root / "dna.synthetic.quickcheck.tsv").read_text(encoding="utf-8"),
                f"id\tbam\texit_code\tstatus\ndna.synthetic\t{bam}\t0\tpass\n",
            )

    def test_idxstats_is_omitted_for_unindexed_rna_contract(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam = self.make_bam(root)
            prefix = root / "rna.synthetic.filtered_cells"

            result = self.run_qc(bam, prefix, False)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((root / "rna.synthetic.filtered_cells.flagstat").exists())
            self.assertTrue((root / "rna.synthetic.filtered_cells.stats").exists())
            self.assertTrue((root / "rna.synthetic.filtered_cells.quickcheck.tsv").exists())
            self.assertFalse((root / "rna.synthetic.filtered_cells.idxstats").exists())

    def test_invalid_bam_records_failed_quickcheck_before_other_qc_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam = self.make_bam(root)
            invalid = root / "invalid.bam"
            invalid.write_text("not a BAM\n", encoding="utf-8")
            prefix = root / "invalid"

            result = self.run_qc(invalid, prefix, False)
            self.assertNotEqual(result.returncode, 0)
            report = (root / "invalid.quickcheck.tsv").read_text(encoding="utf-8")
            self.assertIn("\tfail\n", report)


if __name__ == "__main__":
    unittest.main()
