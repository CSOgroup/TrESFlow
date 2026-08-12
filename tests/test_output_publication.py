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
        rna_workflow = (REPO_ROOT / "subworkflows/local/rna_core/main.nf").read_text(
            encoding="utf-8"
        )
        dna_workflow = (REPO_ROOT / "subworkflows/local/dna_core/main.nf").read_text(
            encoding="utf-8"
        )
        schema = json.loads(
            (REPO_ROOT / "nextflow_schema.json").read_text(encoding="utf-8")
        )
        schema_property = schema["$defs"]["input_output_options"]["properties"][
            "publish_split_fastqs"
        ]

        self.assertIn("publish_split_fastqs = false", nextflow_config)
        self.assertEqual(schema_property["type"], "boolean")
        self.assertFalse(schema_property["default"])
        self.assertNotIn("params.publish_split_fastqs", compression_module)
        self.assertIn("if( params.publish_split_fastqs )", rna_workflow)
        self.assertIn("COMPRESS_RNA_SPLIT_FASTQS(ch_compress_split_input)", rna_workflow)
        self.assertIn("if( params.publish_split_fastqs )", dna_workflow)
        self.assertIn("COMPRESS_DNA_SPLIT_FASTQS(ch_compress_split_input)", dna_workflow)

    def test_report_and_qc_publish_under_tres_stats(self):
        report_module = (
            REPO_ROOT / "modules/local/tres_report_html/main.nf"
        ).read_text(encoding="utf-8")
        samtools_qc_module = (
            REPO_ROOT / "modules/local/samtools_bam_qc/main.nf"
        ).read_text(encoding="utf-8")
        modules_config = (REPO_ROOT / "conf/modules.config").read_text(
            encoding="utf-8"
        )

        self.assertIn('${runtimeOutdir()}/TrES_Stats"', report_module)
        self.assertIn("${runtimeOutdir()}/TrES_Stats/qc/samtools", samtools_qc_module)
        self.assertIn("/TrES_Stats/qc/fastqc", modules_config)
        self.assertIn("/TrES_Stats/qc/samtools", samtools_qc_module)
        self.assertIn("/TrES_Stats/qc/multiqc", modules_config)
        self.assertNotIn("resolvedOutdir')}/qc/", modules_config)
        self.assertNotIn("resolvedOutdir')}/multiqc", modules_config)

    def test_rna_filtered_bam_is_published_directly_without_recompression(self):
        module = (REPO_ROOT / "modules/local/rna_filtered_bam/main.nf").read_text(
            encoding="utf-8"
        )
        workflow = (REPO_ROOT / "subworkflows/local/rna_core/main.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn('${runtimeOutdir()}/rna_align', module)
        self.assertIn('path("${splitName}.filtered_cells.bam")', module)
        self.assertIn("aligned_filtered_bams = RNA_FILTERED_BAM.out.filtered_bam", workflow)
        self.assertNotIn("COMPRESS_RNA_FILTERED_BAM", workflow)
        self.assertFalse((REPO_ROOT / "modules/local/compress_rna_filtered_bam").exists())

    def test_samtools_qc_is_one_combined_process(self):
        workflow = (REPO_ROOT / "workflows/treseq.nf").read_text(encoding="utf-8")
        combined = (REPO_ROOT / "modules/local/samtools_bam_qc/main.nf").read_text(
            encoding="utf-8"
        )

        self.assertIn("include { SAMTOOLS_BAM_QC }", workflow)
        self.assertEqual(workflow.count("SAMTOOLS_BAM_QC(ch_bams_for_samtools_qc)"), 1)
        for obsolete in (
            "SAMTOOLS_FLAGSTAT",
            "SAMTOOLS_STATS",
            "SAMTOOLS_IDXSTATS",
            "SAMTOOLS_QUICKCHECK",
            "SAMTOOLS_QUICKCHECK_REPORT",
        ):
            self.assertNotIn(obsolete, workflow)

        for suffix in ("*.flagstat", "*.stats", "*.idxstats", "*.quickcheck.tsv"):
            self.assertIn(suffix, combined)


if __name__ == "__main__":
    unittest.main()
