import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[1]
MODULES = tuple(
    REPO / relative
    for relative in (
        "modules/local/tag_rna_sb/main.nf",
        "modules/local/tag_rna_umi/main.nf",
        "modules/local/tag_rna_cell_barcode/main.nf",
        "modules/local/split_rna_reads/main.nf",
        "modules/local/fq_to_sam/main.nf",
        "modules/local/tag_dna_sb/main.nf",
        "modules/local/tag_dna_modality/main.nf",
        "modules/local/tag_dna_cell_barcode/main.nf",
        "modules/local/split_dna_reads/main.nf",
    )
)


class Phase2CodonSeqArchitectureTests(unittest.TestCase):
    def test_every_active_process_uses_the_exact_family_environment(self):
        declaration = 'conda "${moduleDir}/../codon_seq/environment.yml"'
        for module in MODULES:
            self.assertIn(declaration, module.read_text(encoding="utf-8"), module)

        environment = (
            REPO / "modules/local/codon_seq/environment.yml"
        ).read_text(encoding="utf-8")
        for requirement in (
            "conda-forge::python=3.12.13=h8ab3286_1_cpython",
            "HCC::codon=0.16.3=haac0b7c_0",
            "HCC::seq=0.11.3=hcb33b08_0",
            "conda-forge::pigz=2.8=h421ea60_2",
            "conda-forge::libxml2=2.11.6=h232c23b_0",
        ):
            self.assertIn(requirement, environment)

    def test_tasks_use_staged_sources_and_normal_path(self):
        forbidden = (
            "runtimeShellExports",
            "runtimeCoreScriptsDir",
            "CODON_BIN",
            "CODON_HOME",
            "PYTHON3_BIN",
            "PIGZ_BIN",
            "${projectDir}",
            "/home/annan",
        )
        for module in MODULES:
            text = module.read_text(encoding="utf-8")
            self.assertIn('path helperScripts, stageAs: "tresflow/bin/*"', text)
            self.assertIn('path codonScripts, stageAs: "tresflow/codon/*"', text)
            self.assertIn('export TMPDIR="\\$PWD/.tmp"', text)
            self.assertIn("python3", text)
            for token in forbidden:
                self.assertNotIn(token, text, f"{token} in {module}")
        for wrapper in (REPO / "bin/tresflow_fastq_utils.py", REPO / "bin/run_fq_to_sam.py"):
            self.assertNotIn("CODON_BIN", wrapper.read_text(encoding="utf-8"))

    def test_production_workflows_supply_every_source_as_path_input(self):
        rna = (REPO / "subworkflows/local/rna_core/main.nf").read_text(encoding="utf-8")
        dna = (REPO / "subworkflows/local/dna_core/main.nf").read_text(encoding="utf-8")
        for name in (
            "run_tag.py",
            "run_tag_lig3.py",
            "tresflow_fastq_utils.py",
            "Tag.codon",
            "Tag_Lig3.codon",
            "utils.codon",
            "Split_ReadsV2.codon",
        ):
            self.assertIn(name, rna + dna)
        for name in ("run_tag_umi.py", "Tag_UMI.codon", "run_fq_to_sam.py", "FqToSAM.codon"):
            self.assertIn(name, rna)
        self.assertNotIn('--script "${coreScriptsDir}', rna + dna)

    def test_global_host_contract_no_longer_requires_codon_or_seq(self):
        runtime = (REPO / "lib/RuntimeSupport.groovy").read_text(encoding="utf-8")
        entry = (REPO / "main.nf").read_text(encoding="utf-8")
        task_exports = (
            REPO / "modules/local/runtime_support/main.nf"
        ).read_text(encoding="utf-8")
        self.assertNotIn("runCodonSeqPreflight", runtime + entry)
        self.assertNotIn("CODON_BIN", runtime + task_exports)
        self.assertNotIn("CODON_HOME", runtime + task_exports)
        self.assertNotIn("codon_home", runtime)
        self.assertNotIn("[name: 'codon', binary: 'codon']", runtime)
        for retained in ("STAR", "samtools", "bwa-mem2", "gatk"):
            self.assertIn(retained, runtime)

    def test_runtime_build_inputs_and_platform_are_immutable(self):
        manifest = json.loads(
            (REPO / "containers/codon-seq/runtime-manifest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(manifest["supported_platforms"], ["linux/amd64"])
        self.assertEqual(
            manifest["publication"]["status"], "local-build-only"
        )
        self.assertEqual(
            manifest["licensing"]["seq_0.11.3"]["upstream_tag_license"],
            "NOASSERTION",
        )
        dockerfile = (REPO / "containers/codon-seq/Dockerfile").read_text(
            encoding="utf-8"
        )
        self.assertIn("FROM --platform=linux/amd64", dockerfile)
        self.assertIn("@sha256:", dockerfile)
        for download in manifest["downloads"]:
            self.assertIn(download["url"], dockerfile)
            self.assertIn(f"--checksum=sha256:{download['sha256']}", dockerfile)
        self.assertIn("codon run -plugin seq -release", dockerfile)

    def test_no_unpublished_image_is_referenced_by_production_processes(self):
        for module in MODULES:
            text = module.read_text(encoding="utf-8")
            self.assertNotIn("container ", text)
            self.assertNotIn("phase2-local", text)
        test_config = (REPO / "tests/phase2_codon_seq.config").read_text(
            encoding="utf-8"
        )
        self.assertIn("tresflow-codon-seq:phase2-local", test_config)
        workflow = (
            REPO / ".github/workflows/codon-seq-runtime.yml"
        ).read_text(encoding="utf-8")
        self.assertNotIn("docker push", workflow)
        self.assertNotIn("push: true", workflow)

    def test_isolated_harness_poisoned_the_legacy_host_prefix(self):
        harness = (REPO / "tests/phase2_codon_seq.nf").read_text(encoding="utf-8")
        self.assertIn("phase2-must-not-be-used", harness)
        for process_name in (
            "TAG_RNA_SAMPLE_BARCODE",
            "TAG_RNA_UMI",
            "TAG_RNA_CELL_BARCODE",
            "SPLIT_RNA_READS",
            "FQ_TO_SAM",
            "TAG_DNA_SAMPLE_BARCODE",
            "TAG_DNA_MODALITY_BARCODE",
            "TAG_DNA_CELL_BARCODE",
            "SPLIT_DNA_READS",
        ):
            self.assertIn(f"{process_name}(", harness)


if __name__ == "__main__":
    unittest.main()
