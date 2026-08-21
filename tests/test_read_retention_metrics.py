import csv
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "bin"))


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SPLIT_RNA = load_module("retention_split_rna", "bin/run_split_reads_rna.py")
SPLIT_DNA = load_module("retention_split_dna", "bin/run_split_reads_dna.py")


def write_fastq(path, comments):
    with path.open("w", encoding="utf-8") as handle:
        for index, comment in enumerate(comments, start=1):
            read_name = f"AV240401:AVT0507:2528453125:1:11104:{index}:3419:UMI{index}"
            handle.write(f"@{read_name} {comment}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")


def read_metrics(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ReadRetentionMetricTests(unittest.TestCase):
    def test_rna_filter_reports_the_existing_nested_predicates(self):
        samtools = shutil.which("samtools")
        if not samtools:
            self.skipTest("samtools is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample = "rna_group"
            source_sam = root / "source.sam"
            sequence = "ACGTACGTACGTACGTACGT"
            quality = "IIIIIIIIIIIIIIIIIIII"
            source_sam.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@SQ\tSN:chrUn\tLN:1000\n"
                "@RG\tID:cell1\tSM:test\n"
                "@RG\tID:other\tSM:test\n"
                f"keep\t99\tchr1\t101\t60\t20M\t=\t151\t70\t{sequence}\t{quality}\tRG:Z:cell1\n"
                f"keep\t147\tchr1\t151\t60\t20M\t=\t101\t-70\t{sequence}\t{quality}\tRG:Z:cell1\n"
                f"other_cell\t99\tchr1\t201\t60\t20M\t=\t251\t70\t{sequence}\t{quality}\tRG:Z:other\n"
                f"other_cell\t147\tchr1\t251\t60\t20M\t=\t201\t-70\t{sequence}\t{quality}\tRG:Z:other\n"
                f"not_paired\t64\tchr1\t301\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:cell1\n"
                f"not_paired\t128\tchr1\t351\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:cell1\n"
                f"noncanonical\t99\tchrUn\t101\t60\t20M\t=\t151\t70\t{sequence}\t{quality}\tRG:Z:cell1\n"
                f"noncanonical\t147\tchrUn\t151\t60\t20M\t=\t101\t-70\t{sequence}\t{quality}\tRG:Z:cell1\n",
                encoding="utf-8",
            )
            aligned_bam = root / f"{sample}.Aligned.sortedByCoord.out.bam"
            subprocess.run(
                [samtools, "sort", "-o", str(aligned_bam), str(source_sam)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run([samtools, "index", str(aligned_bam)], check=True)
            solo = root / "solo" / "filtered"
            solo.mkdir(parents=True)
            (solo / "barcodes.tsv").write_text("cell1\n", encoding="utf-8")
            canonical = root / "canonical.txt"
            canonical.write_text("chr1\n", encoding="utf-8")
            env = os.environ.copy()
            env["SAMTOOLS_BIN"] = samtools
            env["PYTHON3_BIN"] = sys.executable
            subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "scripts/core_runtime/RNA_FILTERED_BAM.sh"),
                    sample,
                    str(root / "solo"),
                    str(aligned_bam),
                    str(canonical),
                    str(root),
                    "1",
                ],
                cwd=root,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )
            values = {
                row["metric"]: int(row["pairs"])
                for row in read_metrics(root / f"{sample}.rna_filter_retention.tsv")
            }
            self.assertEqual(
                values,
                {
                    "star_mapped_primary_pairs": 4,
                    "paired_filter_pairs": 3,
                    "canonical_pairs": 2,
                    "called_cell_pairs": 1,
                },
            )
            filtered_bam = root / f"{sample}.filtered_cells.bam"
            subprocess.run([samtools, "quickcheck", str(filtered_bam)], check=True)
            final_count = subprocess.run(
                [samtools, "view", "--count", "--require-flags", "0x40", str(filtered_bam)],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(final_count.stdout.strip(), "1")

    def test_filter_canonical_combined_validation_preserves_existing_behavior(self):
        samtools = shutil.which("samtools")
        if not samtools:
            self.skipTest("samtools is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_sam = root / "source.sam"
            sequence = "ACGTACGTACGTACGTACGT"
            quality = "IIIIIIIIIIIIIIIIIIII"

            source_sam.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@SQ\tSN:chrUn\tLN:1000\n"
                f"primary\t65\tchr1\t101\t60\t20M\t=\t151\t70\t{sequence}\t{quality}\n"
                f"primary\t129\tchr1\t151\t60\t20M\t=\t101\t-70\t{sequence}\t{quality}\n"
                f"secondary_r1\t321\tchr1\t201\t60\t20M\t=\t251\t70\t{sequence}\t{quality}\n"
                f"supplementary_r1\t2113\tchr1\t301\t60\t20M\t=\t351\t70\t{sequence}\t{quality}\n"
                f"cross_contig_mate\t65\tchr1\t401\t60\t20M\tchrUn\t451\t70\t{sequence}\t{quality}\n"
                f"noncanonical_rname\t0\tchrUn\t501\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\n",
                encoding="utf-8",
            )

            input_bam = root / "input.bam"
            subprocess.run(
                [samtools, "sort", "-o", str(input_bam), str(source_sam)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run([samtools, "index", str(input_bam)], check=True)

            canonical = root / "canonical.txt"
            canonical.write_text("chr1\n", encoding="utf-8")

            output_bam = root / "output.bam"
            summary = root / "validation.tsv"

            env = os.environ.copy()
            env["SAMTOOLS_BIN"] = samtools

            subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "scripts/core_runtime/FilterCanonicalBam.sh"),
                    str(input_bam),
                    str(output_bam),
                    str(canonical),
                    "1",
                    "normal",
                    "--validation-summary",
                    str(summary),
                ],
                cwd=root,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )

            subprocess.run([samtools, "quickcheck", str(output_bam)], check=True)

            # The noncanonical RNAME record must still be excluded.
            retained_rnames = subprocess.run(
                [samtools, "view", str(output_bam)],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.splitlines()
            self.assertTrue(retained_rnames)
            self.assertTrue(
                all(line.split("\t")[2] == "chr1" for line in retained_rnames)
            )

            # A retained chr1 record refers to chrUn as its mate. The helper
            # must therefore preserve chrUn in @SQ exactly as before.
            header = subprocess.run(
                [samtools, "view", "-H", str(output_bam)],
                check=True,
                capture_output=True,
                text=True,
            ).stdout
            self.assertIn("@SQ\tSN:chr1\t", header)
            self.assertIn("@SQ\tSN:chrUn\t", header)

            validation = {}
            with summary.open(encoding="utf-8") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    validation[row["metric"]] = int(row["count"])

            old_count = subprocess.run(
                [
                    samtools,
                    "view",
                    "--count",
                    "--require-flags",
                    "0x40",
                    "--exclude-flags",
                    "0x900",
                    str(output_bam),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertEqual(validation["invalid_alignment_records"], 0)
            self.assertEqual(validation["invalid_mate_reference_records"], 1)
            self.assertEqual(
                validation["primary_r1_records"],
                int(old_count.stdout.strip()),
            )

    def test_filter_canonical_removes_unused_noncanonical_sq_when_mates_are_canonical(self):
        samtools = shutil.which("samtools")
        if not samtools:
            self.skipTest("samtools is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_sam = root / "source.sam"
            sequence = "ACGTACGTACGTACGTACGT"
            quality = "IIIIIIIIIIIIIIIIIIII"

            # chrUn exists in the input dictionary but no retained alignment
            # or mate refers to it. The canonical filter should therefore
            # safely remove chrUn from the output @SQ dictionary.
            source_sam.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@SQ\tSN:chrUn\tLN:1000\n"
                f"canonical_pair\t99\tchr1\t101\t60\t20M\t=\t151\t70\t{sequence}\t{quality}\n"
                f"canonical_pair\t147\tchr1\t151\t60\t20M\t=\t101\t-70\t{sequence}\t{quality}\n"
                f"noncanonical\t0\tchrUn\t501\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\n",
                encoding="utf-8",
            )

            input_bam = root / "input.bam"
            subprocess.run(
                [samtools, "sort", "-o", str(input_bam), str(source_sam)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [samtools, "index", str(input_bam)],
                check=True,
                capture_output=True,
                text=True,
            )

            canonical = root / "canonical.txt"
            canonical.write_text("chr1\n", encoding="utf-8")

            output_bam = root / "output.bam"
            summary = root / "validation.tsv"

            env = os.environ.copy()
            env["SAMTOOLS_BIN"] = samtools

            subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "scripts/core_runtime/FilterCanonicalBam.sh"),
                    str(input_bam),
                    str(output_bam),
                    str(canonical),
                    "1",
                    "normal",
                    "--validation-summary",
                    str(summary),
                ],
                cwd=root,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )

            subprocess.run(
                [samtools, "quickcheck", str(output_bam)],
                check=True,
                capture_output=True,
                text=True,
            )

            header = subprocess.run(
                [samtools, "view", "-H", str(output_bam)],
                check=True,
                capture_output=True,
                text=True,
            ).stdout

            self.assertIn("@SQ\tSN:chr1\t", header)
            self.assertNotIn("@SQ\tSN:chrUn\t", header)

            records = subprocess.run(
                [samtools, "view", str(output_bam)],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.splitlines()

            self.assertEqual(len(records), 2)
            self.assertTrue(
                all(line.split("\t")[2] == "chr1" for line in records)
            )

            validation = {}
            with summary.open(encoding="utf-8") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    validation[row["metric"]] = int(row["count"])

            self.assertEqual(validation["invalid_alignment_records"], 0)
            self.assertEqual(validation["invalid_mate_reference_records"], 0)
            self.assertEqual(validation["primary_r1_records"], 1)

    def test_align_dna_reports_bwa_blacklist_and_proper_pair_populations(self):
        samtools = shutil.which("samtools")
        if not samtools:
            self.skipTest("samtools is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_sam = root / "source.sam"
            source_sam.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "proper\t99\tchr1\t101\t60\t20M\t=\t151\t70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
                "proper\t147\tchr1\t151\t60\t20M\t=\t101\t-70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
                "improper\t65\tchr1\t301\t60\t20M\t=\t351\t70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
                "improper\t129\tchr1\t351\t60\t20M\t=\t301\t-70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
                "blacklisted\t99\tchr1\t501\t60\t20M\t=\t551\t70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
                "blacklisted\t147\tchr1\t551\t60\t20M\t=\t501\t-70\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n",
                encoding="utf-8",
            )
            fake_bwa = root / "fake-bwa-mem2"
            fake_bwa.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "output=''\n"
                "while [[ $# -gt 0 ]]; do\n"
                "  if [[ $1 == -o ]]; then output=$2; shift 2; else shift; fi\n"
                "done\n"
                "cp \"$FAKE_BWA_SAM\" \"$output\"\n",
                encoding="utf-8",
            )
            fake_bwa.chmod(0o755)
            blacklist = root / "blacklist.bed"
            blacklist.write_text("chr1\t490\t600\n", encoding="utf-8")
            rg_header = root / "rg.tsv"
            rg_header.write_text(
                "@RG\tID:sample_H3K27ac\tSM:sample\tLB:test\tPL:ILLUMINA\n",
                encoding="utf-8",
            )
            env = os.environ.copy()
            env.update(
                {
                    "BWA_MEM2_BIN": str(fake_bwa),
                    "FAKE_BWA_SAM": str(source_sam),
                    "SAMTOOLS_BIN": samtools,
                    "ALIGN_DNA_THREADS": "1",
                    "ALIGN_DNA_VIEW_THREADS": "1",
                    "ALIGN_DNA_SORT_THREADS": "1",
                    "ALIGN_DNA_SORT_MEM": "16M",
                }
            )
            subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "scripts/core_runtime/AlignDNA.sh"),
                    "H3K27ac",
                    "sample",
                    str(root / "unused_R1.fastq"),
                    str(root / "unused_R2.fastq"),
                    str(blacklist),
                    str(rg_header),
                    str(root / "unused_reference"),
                    "1000",
                    str(root),
                ],
                cwd=root,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )
            values = {
                row["metric"]: int(row["pairs"])
                for row in read_metrics(root / "sample_H3K27ac.dna_alignment_retention.tsv")
            }
            self.assertEqual(
                values,
                {
                    "bwa_primary_pairs": 3,
                    "post_blacklist_primary_pairs": 2,
                    "post_blacklist_mapped_primary_pairs": 2,
                    "proper_pair_primary_pairs": 1,
                },
            )
            final_count = subprocess.run(
                [samtools, "view", "--count", "--require-flags", "0x40", root / "sample_H3K27ac.bam"],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(final_count.stdout.strip(), "1")

    def run_real_codon_split(self, mode):
        codon = shutil.which("codon")
        if not codon:
            self.skipTest("codon is not installed")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            comments = [
                "CB:Z:AAAACGT\tSB:Z:AAA\tUM:Z:TTTT\tMO:Z:MARKA",
                "CB:Z:NoMatch\tSB:Z:AAA\tUM:Z:TTTT\tMO:Z:MARKA",
            ]
            r1 = root / "r1.fastq"
            r2 = root / "r2.fastq"
            write_fastq(r1, comments)
            write_fastq(r2, comments)
            sb_map = root / "sb.tsv"
            sb_map.write_text("sample\tgroup1\tAAA\n", encoding="utf-8")
            mo_map = root / "mo.tsv"
            mo_map.write_text("sample\tgroup1\tH3K27ac\tMARKA\n", encoding="utf-8")
            output = root / "output"
            output.mkdir()
            command = [
                codon,
                "run",
                "-plugin",
                "seq",
                "-release",
                str(REPO_ROOT / "scripts/core_runtime/Split_ReadsV2.codon"),
                "sample",
                str(output),
                "library",
                mode,
                str(mo_map) if mode == "dna" else "-",
                str(r1),
                str(r2),
                str(sb_map),
            ]
            subprocess.run(command, check=True, capture_output=True, text=True)
            metrics = output / f"sample.{mode}_read_retention.tsv"
            values = {row["metric"]: int(row["pairs"]) for row in read_metrics(metrics)}
            self.assertEqual(values["split_input_pairs"], 2)
            self.assertEqual(values["joint_barcode_accepted_pairs"], 1)
            if mode == "dna":
                header = next(output.glob("SAM_RG_Header_sample_*.tsv")).read_text(
                    encoding="utf-8"
                )
                split_header = next(output.glob("sample_*_R1.fastq")).read_text(
                    encoding="utf-8"
                ).splitlines()[0]
                self.assertIn("@RG\tID:AV240401:AVT0507:2528453125:L1\tSM:sample\tLB:library\tPU:AV240401:AVT0507:2528453125:L1", header)
                self.assertIn("CB:Z:ACGT", split_header)
                self.assertIn("RG:Z:AV240401:AVT0507:2528453125:L1", split_header)

    def test_real_codon_rna_split_emits_retention_metrics(self):
        self.run_real_codon_split("rna")

    def test_real_rna_wrapper_routes_multiple_sb_groups(self):
        codon = shutil.which("codon")
        if not codon:
            self.skipTest("codon is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_id = "multi_group_sample"
            comments = [
                "CB:Z:TAAAAACGT\tSB:Z:TAAAA\tUM:Z:TTTT",
                "CB:Z:GCCCCACGT\tSB:Z:GCCCC\tUM:Z:TTTT",
            ]
            r1 = root / "r1.fastq"
            r2 = root / "r2.fastq"
            write_fastq(r1, comments)
            write_fastq(r2, comments)
            sb_map = root / "rna_sb_group_map.tsv"
            sb_map.write_text(
                "sample\tsb_group\tsb_bc\n"
                f"{sample_id}\talpha\tAAAA\n"
                f"{sample_id}\tbeta\tCCCC\n",
                encoding="utf-8",
            )
            output = root / "output"
            temp_root = root / "runtime_tmp"
            output.mkdir()
            temp_root.mkdir()

            env = os.environ.copy()
            env["CODON_BIN"] = codon
            env["TMPDIR"] = str(temp_root)

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "bin/run_split_reads_rna.py"),
                    "--mode",
                    "real",
                    "--script",
                    str(REPO_ROOT / "scripts/core_runtime/Split_ReadsV2.codon"),
                    "--r1",
                    str(r1),
                    "--r2",
                    str(r2),
                    "--sb-group-map",
                    str(sb_map),
                    "--sample",
                    sample_id,
                    "--library-name",
                    "library",
                    "--output-dir",
                    str(output),
                ],
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )
            rows = read_metrics(output / f"{sample_id}.rna_read_retention.tsv")
            group_counts = {
                row["group"]: int(row["pairs"])
                for row in rows
                if row["metric"] == "routed_group_pairs"
            }
            self.assertEqual(group_counts, {"alpha": 1, "beta": 1})
            self.assertEqual(
                sum(1 for line in (output / f"{sample_id}_alpha_R1.fastq").read_text().splitlines() if line.startswith("@")),
                1,
            )
            self.assertEqual(
                sum(1 for line in (output / f"{sample_id}_beta_R1.fastq").read_text().splitlines() if line.startswith("@")),
                1,
            )
            self.assertEqual(
                sum(1 for line in (output / f"{sample_id}_beta_R2.fastq").read_text().splitlines() if line.startswith("@")),
                1,
            )
            self.assertIn(
                "SB:Z:TAAAA",
                (output / f"{sample_id}_alpha_R1.fastq").read_text(),
            )
            self.assertIn(
                "SB:Z:GCCCC",
                (output / f"{sample_id}_beta_R1.fastq").read_text(),
            )

    def test_real_codon_dna_split_emits_retention_metrics(self):
        self.run_real_codon_split("dna")

    def run_real_wrapper_split(self, mode):
        codon = shutil.which("codon")
        if not codon:
            self.skipTest("codon is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            temp_root = root / "runtime_tmp"
            output = root / "nextflow_workdir"
            temp_root.mkdir()
            output.mkdir()

            comments = [
                "CB:Z:AAAACGT\tSB:Z:AAA\tUM:Z:TTTT\tMO:Z:MARKA",
                "CB:Z:NoMatch\tSB:Z:AAA\tUM:Z:TTTT\tMO:Z:MARKA",
            ]
            r1 = root / "r1.fastq"
            r2 = root / "r2.fastq"
            write_fastq(r1, comments)
            write_fastq(r2, comments)
            sb_map = root / "sb.tsv"
            sb_map.write_text("sample\tgroup1\tAAA\n", encoding="utf-8")
            mo_map = root / "mo.tsv"
            mo_map.write_text("sample\tgroup1\tH3K27ac\tMARKA\n", encoding="utf-8")

            wrapper = REPO_ROOT / f"bin/run_split_reads_{mode}.py"
            command = [
                sys.executable,
                str(wrapper),
                "--mode",
                "real",
                "--script",
                str(REPO_ROOT / "scripts/core_runtime/Split_ReadsV2.codon"),
                "--r1",
                str(r1),
                "--r2",
                str(r2),
            ]
            if mode == "dna":
                command.extend(["--mo-map", str(mo_map)])
            command.extend(
                [
                    "--sb-group-map",
                    str(sb_map),
                    "--sample",
                    "sample",
                    "--library-name",
                    "library",
                    "--output-dir",
                    str(output),
                ]
            )

            env = os.environ.copy()
            env["CODON_BIN"] = codon
            env["TMPDIR"] = str(temp_root)
            result = subprocess.run(
                command,
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )

            metrics = output / f"sample.{mode}_read_retention.tsv"
            self.assertTrue(metrics.is_file())
            values = {row["metric"]: int(row["pairs"]) for row in read_metrics(metrics)}
            self.assertEqual(values["split_input_pairs"], 2)
            self.assertEqual(values["joint_barcode_accepted_pairs"], 1)
            self.assertIn(f"Finished split output move | {metrics}", result.stderr)
            self.assertEqual(list(temp_root.iterdir()), [])

            self.assertTrue(list(output.glob("sample_*_R1.fastq")))
            self.assertTrue(list(output.glob("sample_*_R2.fastq")))
            self.assertTrue(list(output.glob("SAM_RG_Header_sample_*.tsv")))

    def test_real_rna_wrapper_preserves_retention_metrics_before_temp_cleanup(self):
        self.run_real_wrapper_split("rna")

    def test_real_dna_wrapper_preserves_retention_metrics_before_temp_cleanup(self):
        self.run_real_wrapper_split("dna")

    def test_mock_rna_split_reports_input_joint_and_routed_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            comments = [
                "CB:Z:AAAACGT\tSB:Z:AAAA\tUM:Z:TTTT",
                "CB:Z:NoMatch\tSB:Z:AAAA\tUM:Z:TTTT",
                "CB:Z:CCCCACGT\tSB:Z:CCCC\tUM:Z:TTTT",
            ]
            r1 = root / "r1.fastq"
            r2 = root / "r2.fastq"
            write_fastq(r1, comments)
            write_fastq(r2, comments)
            sb_map = root / "sb.tsv"
            sb_map.write_text("sample\tgroup1\tAAAA\nsample\tgroup2\tCCCC\n", encoding="utf-8")

            SPLIT_RNA.mock_split(
                SimpleNamespace(
                    r1=r1,
                    r2=r2,
                    sb_group_map=sb_map,
                    sample="sample",
                    library_name="library",
                    output_dir=root,
                )
            )

            rows = read_metrics(root / "sample.rna_read_retention.tsv")
            values = {(row["metric"], row["group"]): int(row["pairs"]) for row in rows}
            self.assertEqual(values[("split_input_pairs", "__all__")], 3)
            self.assertEqual(values[("joint_barcode_accepted_pairs", "__all__")], 2)
            self.assertEqual(values[("routed_group_pairs", "group1")], 1)
            self.assertEqual(values[("routed_group_pairs", "group2")], 1)

    def test_mock_dna_split_reports_mark_branch_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            comments = [
                "CB:Z:AAAACGT\tSB:Z:AAA\tMO:Z:MARKA",
                "CB:Z:AAAACGT\tSB:Z:AAA\tMO:Z:MARKB",
                "CB:Z:NoMatch\tSB:Z:AAA\tMO:Z:MARKA",
            ]
            r1 = root / "r1.fastq"
            r2 = root / "r2.fastq"
            write_fastq(r1, comments)
            write_fastq(r2, comments)
            sb_map = root / "sb.tsv"
            sb_map.write_text("sample\tgroup1\tAAA\n", encoding="utf-8")
            mo_map = root / "mo.tsv"
            mo_map.write_text(
                "sample\tgroup1\tH3K27ac\tMARKA\n"
                "sample\tgroup1\tH3K27me3\tMARKB\n",
                encoding="utf-8",
            )

            SPLIT_DNA.mock_split(
                SimpleNamespace(
                    r1=r1,
                    r2=r2,
                    sb_group_map=sb_map,
                    mo_map=mo_map,
                    sample="sample",
                    library_name="library",
                    output_dir=root,
                )
            )

            rows = read_metrics(root / "sample.dna_read_retention.tsv")
            branch = {
                row["branch"]: int(row["pairs"])
                for row in rows
                if row["metric"] == "routed_branch_pairs"
            }
            self.assertEqual(branch, {"H3K27ac": 1, "H3K27me3": 1})

    def test_rna_bam_predicate_audit_is_strictly_nested(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            canonical = root / "canonical.txt"
            canonical.write_text("chr1\n", encoding="utf-8")
            cells = root / "barcodes.tsv"
            cells.write_text("cell1\n", encoding="utf-8")
            output = root / "metrics.tsv"
            sam = "\n".join(
                [
                    "p1\t99\tchr1\t1\t60\t20M\t=\t100\t0\tACGT\tIIII\tRG:Z:cell1",
                    "p2\t65\tchr1\t2\t60\t20M\t*\t0\t0\tACGT\tIIII\tRG:Z:other",
                    "p3\t65\tchrUn\t3\t60\t20M\t*\t0\t0\tACGT\tIIII\tRG:Z:cell1",
                    "p4\t64\tchr1\t4\t60\t20M\t*\t0\t0\tACGT\tIIII\tRG:Z:cell1",
                    "p5\t69\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:cell1",
                ]
            ) + "\n"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts/core_runtime/SummarizeRnaRetention.py"),
                    "--split-id",
                    "sample_group",
                    "--canonical-contigs",
                    str(canonical),
                    "--called-barcodes",
                    str(cells),
                    "--output",
                    str(output),
                ],
                input=sam,
                text=True,
                check=True,
            )

            values = {row["metric"]: int(row["pairs"]) for row in read_metrics(output)}
            self.assertEqual(
                values,
                {
                    "star_mapped_primary_pairs": 4,
                    "paired_filter_pairs": 3,
                    "canonical_pairs": 2,
                    "called_cell_pairs": 1,
                },
            )


if __name__ == "__main__":
    unittest.main()
