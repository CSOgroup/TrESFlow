import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "bin"))

import tresflow_fastq_utils as FASTQ_UTILS


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REPORT = load_module("aviti_report", "bin/render_tres_report.py")
AVITI_REGEX = r"^(?:[^:]+:){4}([0-9]+):([0-9]+):([0-9]+):[^:]+$"
CELL_BARCODE = "ACGTACGTTGCATGCAGATCGATC"


def read_picard_metrics(path):
    header = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line or raw_line.startswith("#"):
            continue
        fields = raw_line.split("\t")
        if header is None:
            if "READ_PAIR_DUPLICATES" in fields:
                header = fields
            continue
        return dict(zip(header, fields))
    raise AssertionError(f"No Picard metrics row found in {path}")


class AvitiReadGroupTests(unittest.TestCase):
    def test_aviti_qname_parsing(self):
        qname = "AV240401:AVT0507:2528453125:1:11104:5031:3419:ACGTACGT"
        self.assertEqual(FASTQ_UTILS.parse_aviti_qname(qname), (1, 11104, 5031, 3419))

    def test_dna_read_groups_are_lane_specific_while_cb_is_unchanged(self):
        comment = (
            f"CB:Z:AAA{CELL_BARCODE}\tRG:Z:AAA{CELL_BARCODE}"
            "\tMO:Z:AGGCTATA\tSB:Z:AAA"
        )
        lane1 = FASTQ_UTILS.canonicalize_dna_fastq_comment(
            "sample", "group", "AV240401:AVT0507:FC:1:11104:5031:3419:UMI1", comment
        )
        lane2 = FASTQ_UTILS.canonicalize_dna_fastq_comment(
            "sample", "group", "AV240401:AVT0507:FC:2:11104:5031:3419:UMI2", comment
        )

        self.assertEqual(FASTQ_UTILS.find_tag_value(lane1, "CB"), CELL_BARCODE)
        self.assertEqual(FASTQ_UTILS.find_tag_value(lane2, "CB"), CELL_BARCODE)
        self.assertEqual(FASTQ_UTILS.find_tag_value(lane1, "RG"), f"{CELL_BARCODE}_L1")
        self.assertEqual(FASTQ_UTILS.find_tag_value(lane2, "RG"), f"{CELL_BARCODE}_L2")

    def test_lane_read_group_headers_share_one_library(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            header = Path(tmpdir) / "rg.tsv"
            FASTQ_UTILS.write_rg_header(
                header,
                "sample",
                "logical_library",
                {f"{CELL_BARCODE}_L1", f"{CELL_BARCODE}_L2"},
            )
            rows = [dict(field.split(":", 1) for field in line.split("\t")[1:])
                    for line in header.read_text(encoding="utf-8").splitlines()]

        self.assertEqual({row["ID"] for row in rows}, {f"{CELL_BARCODE}_L1", f"{CELL_BARCODE}_L2"})
        self.assertEqual({row["LB"] for row in rows}, {"logical_library"})

    def test_parameter_and_markduplicates_configuration(self):
        nextflow_config = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")
        schema = json.loads((REPO_ROOT / "nextflow_schema.json").read_text(encoding="utf-8"))
        module_config = (REPO_ROOT / "conf/modules.config").read_text(encoding="utf-8")
        parameter = schema["$defs"]["execution_options"]["properties"][
            "aviti_optical_duplicate_distance"
        ]

        self.assertIn("aviti_optical_duplicate_distance = 10", nextflow_config)
        self.assertEqual(parameter["default"], 10)
        self.assertEqual(parameter["minimum"], 0)
        self.assertIn("--BARCODE_TAG CB", module_config)
        self.assertIn("--REMOVE_DUPLICATES false", module_config)
        self.assertIn("--READ_NAME_REGEX", module_config)
        self.assertIn(AVITI_REGEX[:-1], module_config)
        self.assertIn("--OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.aviti_optical_duplicate_distance}", module_config)

    def test_report_parser_exposes_complexity_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metrics = Path(tmpdir) / "sample.DuplicateMetrics.txt"
            metrics.write_text(
                "## mock\n"
                "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
                "UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
                "READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n"
                "lib\t0\t100\t0\t20\t15\t0.2\t345\n",
                encoding="utf-8",
            )
            parsed = REPORT.read_duplicate_metrics(metrics)

        self.assertEqual(parsed["read_pairs_examined"], 100)
        self.assertEqual(parsed["read_pair_duplicates"], 20)
        self.assertEqual(parsed["read_pair_optical_duplicates"], 15)
        self.assertEqual(parsed["estimated_library_size"], 345)

    def test_duplicate_removal_remains_flag_0x400_only(self):
        script = (REPO_ROOT / "scripts/core_runtime/SplitDuplicatesDNA.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("--exclude-flags 0x400", script)
        self.assertNotIn("OPTICAL", script)


class PicardAvitiIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gatk = shutil.which("gatk")
        cls.samtools = shutil.which("samtools")

    def make_synthetic_bam(self, root):
        sam = root / "input.sam"
        bam = root / "input.bam"
        sequence = "A" * 50
        quality = "I" * 50
        lines = [
            "@HD\tVN:1.6\tSO:unsorted",
            "@SQ\tSN:chr1\tLN:100000",
            f"@RG\tID:{CELL_BARCODE}_L1\tSM:sample\tLB:logical_library\tPL:ELEMENT",
            f"@RG\tID:{CELL_BARCODE}_L2\tSM:sample\tLB:logical_library\tPL:ELEMENT",
        ]

        def add_pair(qname, lane, first_position, mate_position):
            rg = f"{CELL_BARCODE}_L{lane}"
            template_length = mate_position + 49 - first_position + 1
            tags = f"RG:Z:{rg}\tCB:Z:{CELL_BARCODE}"
            lines.append(
                f"{qname}\t99\tchr1\t{first_position}\t60\t50M\t=\t{mate_position}\t"
                f"{template_length}\t{sequence}\t{quality}\t{tags}"
            )
            lines.append(
                f"{qname}\t147\tchr1\t{mate_position}\t60\t50M\t=\t{first_position}\t"
                f"-{template_length}\t{sequence}\t{quality}\t{tags}"
            )

        # One genomic duplicate family: two lane-1 members are spatial neighbors;
        # the lane-2 member has colliding tile/x/y but is physically independent.
        add_pair("AV240401:AVT0507:FC:1:11104:100:100:UMI1", 1, 101, 201)
        add_pair("AV240401:AVT0507:FC:1:11104:105:106:UMI2", 1, 101, 201)
        add_pair("AV240401:AVT0507:FC:2:11104:105:106:UMI3", 2, 101, 201)

        # Unique pairs keep the library-size estimate away from tiny-sample rounding.
        for index in range(100):
            first = 1000 + (index * 300)
            add_pair(
                f"AV240401:AVT0507:FC:1:{12000 + index}:1000:1000:U{index}",
                1,
                first,
                first + 100,
            )

        sam.write_text("\n".join(lines) + "\n", encoding="utf-8")
        subprocess.run(
            [self.samtools, "sort", "-o", str(bam), str(sam)],
            check=True,
            capture_output=True,
            text=True,
        )
        return bam

    def run_markduplicates(self, root, input_bam, cutoff):
        output = root / f"marked_{cutoff}.bam"
        metrics = root / f"marked_{cutoff}.metrics"
        command = [
            self.gatk,
            "--java-options",
            "-Xmx1g -XX:-UsePerfData",
            "MarkDuplicates",
            "--INPUT",
            str(input_bam),
            "--OUTPUT",
            str(output),
            "--METRICS_FILE",
            str(metrics),
            "--REMOVE_DUPLICATES",
            "false",
            "--BARCODE_TAG",
            "CB",
            "--READ_NAME_REGEX",
            AVITI_REGEX,
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE",
            str(cutoff),
            "--CREATE_INDEX",
            "false",
            "--VALIDATION_STRINGENCY",
            "SILENT",
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)
        return output, read_picard_metrics(metrics)

    def test_picard_groups_genomic_duplicates_cross_lane_but_optical_duplicates_within_lane(self):
        if not self.gatk or not self.samtools:
            self.skipTest("gatk and samtools are required for the targeted integration test")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_bam = self.make_synthetic_bam(root)
            output10, metrics10 = self.run_markduplicates(root, input_bam, 10)
            _, metrics0 = self.run_markduplicates(root, input_bam, 0)
            duplicate_pairs = subprocess.run(
                [
                    self.samtools,
                    "view",
                    "--count",
                    "--require-flags",
                    "0x440",
                    str(output10),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

        # Shared LB + CB keeps all three coordinate-identical molecules in one
        # genomic family, so two read pairs are marked duplicate across lanes.
        self.assertEqual(int(metrics10["READ_PAIR_DUPLICATES"]), 2)
        self.assertEqual(int(duplicate_pairs.stdout.strip()), 2)
        # RG separates lanes physically: only the two lane-1 members cluster.
        self.assertEqual(int(metrics10["READ_PAIR_OPTICAL_DUPLICATES"]), 1)
        self.assertEqual(int(metrics0["READ_PAIR_OPTICAL_DUPLICATES"]), 0)
        self.assertGreater(
            int(metrics10["ESTIMATED_LIBRARY_SIZE"]),
            int(metrics0["ESTIMATED_LIBRARY_SIZE"]),
        )


if __name__ == "__main__":
    unittest.main()
