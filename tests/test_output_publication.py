import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class OutputPublicationTests(unittest.TestCase):
    def test_split_fastq_publication_defaults_off_and_gates_compression(self):
        nextflow_config = (REPO_ROOT / "nextflow.config").read_text(encoding="utf-8")
        compression_module = (
            REPO_ROOT / "modules/local/compress_split_fastqs/main.nf"
        ).read_text(encoding="utf-8")
        schema = json.loads(
            (REPO_ROOT / "nextflow_schema.json").read_text(encoding="utf-8")
        )
        schema_property = schema["$defs"]["input_output_options"]["properties"][
            "publish_split_fastqs"
        ]

        self.assertIn("publish_split_fastqs = false", nextflow_config)
        self.assertEqual(schema_property["type"], "boolean")
        self.assertFalse(schema_property["default"])
        self.assertIn("when:\n    params.publish_split_fastqs", compression_module)

    def test_report_and_qc_publish_under_tres_stats(self):
        report_module = (
            REPO_ROOT / "modules/local/tres_report_html/main.nf"
        ).read_text(encoding="utf-8")
        quickcheck_module = (
            REPO_ROOT / "modules/local/samtools_quickcheck_report/main.nf"
        ).read_text(encoding="utf-8")
        modules_config = (REPO_ROOT / "conf/modules.config").read_text(
            encoding="utf-8"
        )

        self.assertIn('${runtimeOutdir()}/TrES_Stats"', report_module)
        self.assertIn("${runtimeOutdir()}/TrES_Stats/qc/samtools", quickcheck_module)
        self.assertIn("/TrES_Stats/qc/fastqc", modules_config)
        self.assertIn("/TrES_Stats/qc/samtools", modules_config)
        self.assertIn("/TrES_Stats/qc/multiqc", modules_config)
        self.assertNotIn("resolvedOutdir')}/qc/", modules_config)
        self.assertNotIn("resolvedOutdir')}/multiqc", modules_config)


if __name__ == "__main__":
    unittest.main()
