import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SPLIT_SCRIPT = REPO_ROOT / "scripts/core_runtime/SplitDuplicatesDNA.sh"
CANONICAL_FILTER = REPO_ROOT / "scripts/core_runtime/FilterCanonicalBam.sh"


class DnaNoDupOptimizationTests(unittest.TestCase):
    def setUp(self):
        self.samtools = shutil.which("samtools")
        if not self.samtools:
            self.skipTest("samtools is not installed")

    def make_bam(self, root, name, records):
        sam = root / f"{name}.sam"
        header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "@SQ\tSN:chrX\tLN:1000\n"
            "@SQ\tSN:chrY\tLN:1000\n"
            "@SQ\tSN:chrM\tLN:1000\n"
        )
        sam.write_text(header + "".join(records), encoding="utf-8")
        bam = root / f"{name}.bam"
        subprocess.run(
            [self.samtools, "view", "--bam", "--output", str(bam), str(sam)],
            check=True,
        )
        subprocess.run([self.samtools, "index", str(bam)], check=True)
        return bam

    @staticmethod
    def mapped_pair(name, contig, position, duplicate=False):
        duplicate_bit = 0x400 if duplicate else 0
        first_flag = 99 | duplicate_bit
        second_flag = 147 | duplicate_bit
        return [
            f"{name}\t{first_flag}\t{contig}\t{position}\t60\t20M\t=\t{position + 50}\t70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n",
            f"{name}\t{second_flag}\t{contig}\t{position + 50}\t60\t20M\t=\t{position}\t-70\tTGCATGCATGCATGCATGCA\tIIIIIIIIIIIIIIIIIIII\n",
        ]

    def run_split(self, root, input_bam, split_name):
        output_bam = root / f"{split_name}_NoDup.bam"
        output_bai = root / f"{split_name}_NoDup.bam.bai"
        mapped_reads = root / f"{split_name}.nodup_mapped_reads.txt"
        warning = root / f"{split_name}.zero_mapped_nodup_bam.tsv"
        env = os.environ.copy()
        env["SAMTOOLS_BIN"] = self.samtools
        subprocess.run(
            [
                "bash",
                str(SPLIT_SCRIPT),
                str(input_bam),
                str(output_bam),
                str(output_bai),
                str(mapped_reads),
                str(warning),
                "2",
                "sample",
                "group",
                "mark",
                split_name,
            ],
            check=True,
            env=env,
            capture_output=True,
            text=True,
        )
        return output_bam, output_bai, mapped_reads, warning

    def samtools_view(self, bam):
        return subprocess.run(
            [self.samtools, "view", str(bam)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout

    def test_duplicate_filter_is_equivalent_to_former_canonical_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            records = []
            records += self.mapped_pair("keep_chr1", "chr1", 101)
            records += self.mapped_pair("drop_duplicate", "chr1", 301, duplicate=True)
            records += self.mapped_pair("keep_chrX", "chrX", 101)
            records += self.mapped_pair("keep_chrY", "chrY", 101)
            records += self.mapped_pair("keep_chrM", "chrM", 101)
            markeddup = self.make_bam(root, "canonical_markeddup", records)

            output_bam, output_bai, mapped_reads, warning = self.run_split(
                root, markeddup, "sample_group_mark"
            )

            allowlist = root / "canonical.txt"
            allowlist.write_text("chr1\nchrX\nchrY\nchrM\n", encoding="utf-8")
            baseline = root / "former_NoDup.bam"
            env = os.environ.copy()
            env.update(
                {
                    "SAMTOOLS_BIN": self.samtools,
                    "TMPDIR": str(root),
                }
            )
            subprocess.run(
                [
                    "bash",
                    str(CANONICAL_FILTER),
                    str(markeddup),
                    str(baseline),
                    str(allowlist),
                    "2",
                    "normal",
                    "--exclude-flags",
                    "0x400",
                ],
                check=True,
                env=env,
            )

            self.assertEqual(self.samtools_view(output_bam), self.samtools_view(baseline))
            self.assertEqual(mapped_reads.read_text(encoding="utf-8"), "8\n")
            self.assertFalse(warning.exists())
            subprocess.run([self.samtools, "quickcheck", str(output_bam)], check=True)
            self.assertTrue(output_bai.exists())

            output_records = self.samtools_view(output_bam).splitlines()
            self.assertEqual(
                {record.split("\t")[2] for record in output_records},
                {"chr1", "chrX", "chrY", "chrM"},
            )
            self.assertFalse(
                any(int(record.split("\t")[1]) & 0x400 for record in output_records)
            )
            for contig in ("chr1", "chrX", "chrY", "chrM"):
                subprocess.run(
                    [self.samtools, "view", str(output_bam), contig],
                    check=True,
                    capture_output=True,
                    text=True,
                )

    def test_zero_mapped_bam_is_indexed_and_writes_historical_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            records = [
                "unmapped\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
                "unmapped\t141\t*\t0\t0\t*\t*\t0\t0\tTGCA\tIIII\n",
            ]
            markeddup = self.make_bam(root, "zero_mapped_markeddup", records)
            output_bam, output_bai, mapped_reads, warning = self.run_split(
                root, markeddup, "sample_group_mark"
            )

            self.assertEqual(mapped_reads.read_text(encoding="utf-8"), "0\n")
            self.assertTrue(output_bai.exists())
            subprocess.run([self.samtools, "quickcheck", str(output_bam)], check=True)
            expected = (
                "sample\tgroup\tmark\tbam\tmapped_reads\tskipped_output\n"
                f"sample\tgroup\tmark\t{output_bam.resolve()}\t0\t"
                "sample_group_mark_NoDup.bw\n"
            )
            self.assertEqual(warning.read_text(encoding="utf-8"), expected)

    def test_workflow_has_no_bam_copy_gate_or_coverage_normalizer(self):
        workflow = (REPO_ROOT / "subworkflows/local/dna_core/main.nf").read_text(
            encoding="utf-8"
        )
        split_module = (
            REPO_ROOT / "modules/local/split_duplicates_dna/main.nf"
        ).read_text(encoding="utf-8")
        split_script = SPLIT_SCRIPT.read_text(encoding="utf-8")
        coverage_module = (
            REPO_ROOT / "modules/nf-core/deeptools/bamcoverage/main.nf"
        ).read_text(encoding="utf-8")

        self.assertNotIn("CHECK_DNA_NODUP_BAM", workflow)
        self.assertNotIn("NORMALIZE_DNA_BAMCOVERAGE", workflow)
        self.assertIn("SPLIT_DUPLICATES_DNA.out.mapped_reads", workflow)
        self.assertIn("coverage_warnings = SPLIT_DUPLICATES_DNA.out.warnings", workflow)
        self.assertIn("coverage_bigwigs = ch_coverage_bigwigs", workflow)
        self.assertNotIn("canonicalChromosomes", split_module)
        self.assertNotIn("FilterCanonicalBam.sh", split_module)
        self.assertNotIn("cp ", split_script)
        self.assertIn("--exclude-flags 0x400", split_script)
        self.assertIn('path("*.bw")', coverage_module)
        self.assertNotIn("bigWig", coverage_module)
        self.assertFalse((REPO_ROOT / "modules/local/check_dna_nodup_bam").exists())
        self.assertFalse((REPO_ROOT / "modules/local/normalize_dna_bamcoverage").exists())


if __name__ == "__main__":
    unittest.main()
