import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[1]
TRIM_MODULES = (
    REPO / "modules/local/trim_rna_fastqs/main.nf",
    REPO / "modules/local/trim_dna_fastqs/main.nf",
)
FILTER_MODULE = REPO / "modules/local/dual_tag_artifact_filter/main.nf"
COMPRESS_MODULE = REPO / "modules/local/compress_split_fastqs/main.nf"
TRIM_CONTAINER = (
    "quay.io/biocontainers/trim-galore@"
    "sha256:a02bb87b8ce02d86efd0ffd65e2cce1559b52689faab42faad1df145657390cf"
)
CUTADAPT_CONTAINER = (
    "quay.io/biocontainers/cutadapt@"
    "sha256:2049f305574854edb189ccad7038fda4801ef16458bcde0239383d42d4a3f83a"
)


class Phase3FastqPreprocessingArchitectureTests(unittest.TestCase):
    def test_trim_processes_use_exact_environment_container_and_staged_helper(self):
        for module in TRIM_MODULES:
            text = module.read_text(encoding="utf-8")
            self.assertIn("environment-trim.yml", text)
            self.assertIn(TRIM_CONTAINER, text)
            self.assertIn("path helperScript, stageAs: 'tresflow/bin/run_trim_galore.py'", text)
            self.assertIn('python3 "tresflow/bin/run_trim_galore.py"', text)
            self.assertIn('export TMPDIR="\\$PWD/.tmp"', text)
            self.assertNotIn("codon_wrapper", text)
            for forbidden in (
                "runtimeShellExports",
                "PYTHON3_BIN",
                "TRIM_GALORE_BIN",
                "${projectDir}",
                "/home/annan",
            ):
                self.assertNotIn(forbidden, text, f"{forbidden} in {module}")

    def test_filter_and_compression_use_path_tools(self):
        filter_text = FILTER_MODULE.read_text(encoding="utf-8")
        self.assertIn("environment-cutadapt.yml", filter_text)
        self.assertIn(CUTADAPT_CONTAINER, filter_text)
        self.assertIn(
            "path helperScript, stageAs: 'tresflow/bin/run_dual_tag_artifact_filter.py'",
            filter_text,
        )
        self.assertIn('python3 "tresflow/bin/run_dual_tag_artifact_filter.py"', filter_text)
        self.assertNotIn("--cutadapt-bin", filter_text)
        self.assertIn("cutadapt --version", filter_text)

        compression = COMPRESS_MODULE.read_text(encoding="utf-8")
        self.assertIn("environment-pigz.yml", compression)
        self.assertIn(CUTADAPT_CONTAINER, compression)
        self.assertIn('pigz -c -p "${task.cpus}"', compression)
        self.assertNotIn("PIGZ_BIN", compression)
        for text in (filter_text, compression):
            self.assertIn('export TMPDIR="\\$PWD/.tmp"', text)
            self.assertNotIn("runtimeShellExports", text)
            self.assertNotIn("${projectDir}", text)

    def test_exact_conda_and_container_manifest(self):
        runtime_root = REPO / "modules/local/fastq_preprocessing"
        expected = {
            "environment-trim.yml": (
                "conda-forge::python=3.12.13=hd63d673_0_cpython",
                "bioconda::trim-galore=0.6.11=hdfd78af_0",
                "bioconda::cutadapt=5.2=py312h0fa9677_0",
            ),
            "environment-cutadapt.yml": (
                "conda-forge::python=3.12.13=hd63d673_0_cpython",
                "bioconda::cutadapt=5.2=py312h0fa9677_0",
            ),
            "environment-pigz.yml": ("conda-forge::pigz=2.8=h421ea60_2",),
        }
        for filename, requirements in expected.items():
            environment = (runtime_root / filename).read_text(encoding="utf-8")
            for requirement in requirements:
                self.assertIn(requirement, environment)

        manifest = json.loads((runtime_root / "runtime-manifest.json").read_text())
        self.assertEqual(manifest["containers"]["trim"]["reference"], TRIM_CONTAINER)
        self.assertEqual(
            manifest["containers"]["cutadapt_and_compression"]["reference"],
            CUTADAPT_CONTAINER,
        )
        self.assertEqual(
            {package["name"] for package in manifest["conda_direct_packages"]},
            {"python", "trim-galore", "cutadapt", "pigz"},
        )
        for package in manifest["conda_direct_packages"]:
            self.assertRegex(package["sha256"], r"^[0-9a-f]{64}$")

    def test_production_workflows_stage_helpers_and_cover_every_active_process(self):
        rna = (REPO / "subworkflows/local/rna_core/main.nf").read_text(encoding="utf-8")
        dna = (REPO / "subworkflows/local/dna_core/main.nf").read_text(encoding="utf-8")
        self.assertIn("TRIM_RNA_FASTQS(TAG_RNA_CELL_BARCODE.out.tagged, trimHelperScript)", rna)
        self.assertIn("TRIM_DNA_FASTQS(TAG_DNA_CELL_BARCODE.out.tagged, trimHelperScript)", dna)
        self.assertIn("dualTagFilterHelperScript", dna)
        for process_name in (
            "TRIM_RNA_FASTQS",
            "TRIM_DNA_FASTQS",
            "DUAL_TAG_ARTIFACT_FILTER",
            "COMPRESS_RNA_SPLIT_FASTQS",
            "COMPRESS_DNA_SPLIT_FASTQS",
        ):
            self.assertIn(process_name, rna + dna)

    def test_global_host_contract_no_longer_requires_preprocessing_tools(self):
        runtime = (REPO / "lib/RuntimeSupport.groovy").read_text(encoding="utf-8")
        task_exports = (REPO / "modules/local/runtime_support/main.nf").read_text(
            encoding="utf-8"
        )
        for token in (
            "CUTADAPT_BIN",
            "TRIM_GALORE_BIN",
            "PIGZ_BIN",
            "[name: 'cutadapt', binary: 'cutadapt']",
            "[name: 'trim_galore', binary: 'trim_galore']",
            "[name: 'pigz', binary: 'pigz']",
        ):
            self.assertNotIn(token, runtime + task_exports)
        for retained in ("python3", "samtools", "bwa-mem2", "bamCoverage", "gatk"):
            self.assertIn(retained, runtime)

    def test_wrappers_do_not_read_absolute_host_binary_variables(self):
        trim_wrapper = (REPO / "bin/run_trim_galore.py").read_text(encoding="utf-8")
        filter_wrapper = (REPO / "bin/run_dual_tag_artifact_filter.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn('os.environ["TRIM_GALORE_BIN"]', trim_wrapper)
        self.assertIn('shutil.which("trim_galore")', trim_wrapper)
        self.assertIn('shutil.which("cutadapt")', filter_wrapper)

    def test_fixture_manifest_hashes_and_behavior_are_explicit(self):
        fixture_root = REPO / "tests/fixtures/fastq_preprocessing"
        manifest = json.loads((fixture_root / "manifest.json").read_text())
        self.assertEqual(manifest["asset_type"], "synthetic-test-only")
        self.assertEqual(manifest["license"], "MIT")
        self.assertEqual(manifest["expected_behavior"]["dual_tag_rejected_pairs"], 1)
        for filename, expected_hash in manifest["files"].items():
            observed = hashlib.sha256((fixture_root / filename).read_bytes()).hexdigest()
            self.assertEqual(observed, expected_hash, filename)

    def test_apptainer_profile_uses_the_same_digest_pinned_process_images(self):
        config = (REPO / "tests/phase3_fastq_preprocessing.config").read_text(
            encoding="utf-8"
        )
        self.assertIn("phase3_fastq_apptainer_static", config)
        self.assertIn("apptainer.enabled = true", config)
        for module in (*TRIM_MODULES, FILTER_MODULE, COMPRESS_MODULE):
            self.assertIn("@sha256:", module.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
