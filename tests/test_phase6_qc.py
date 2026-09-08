import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[1]
FASTQC = REPO / "modules/nf-core/fastqc/main.nf"
SAMTOOLS_QC = REPO / "modules/local/samtools_bam_qc/main.nf"
MULTIQC = REPO / "modules/nf-core/multiqc/main.nf"
WORKFLOW = REPO / "workflows/treseq.nf"

CONTAINERS = {
    FASTQC: (
        "quay.io/biocontainers/fastqc@"
        "sha256:e194048df39c3145d9b4e0a14f4da20b59d59250465b6f2a9cb698445fd45900"
    ),
    SAMTOOLS_QC: (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    MULTIQC: (
        "quay.io/biocontainers/multiqc@"
        "sha256:dfd9fde2c48b896b884e79a71ddc16c72c97a0ee5c5c8e45aaba50f55d07d263"
    ),
}


class Phase6QcArchitectureTests(unittest.TestCase):
    def test_active_qc_processes_use_exact_environments_and_containers(self):
        environment_checks = {
            FASTQC: (
                "environment.yml",
                "bioconda::fastqc=0.12.1=hdfd78af_0",
            ),
            SAMTOOLS_QC: (
                "../dna_processing/environment-samtools.yml",
                "bioconda::samtools=1.23.1=ha83d96e_0",
            ),
            MULTIQC: (
                "environment.yml",
                "bioconda::multiqc=1.33=pyhdfd78af_0",
            ),
        }
        environment_paths = {
            FASTQC: FASTQC.parent / "environment.yml",
            SAMTOOLS_QC: REPO / "modules/local/dna_processing/environment-samtools.yml",
            MULTIQC: MULTIQC.parent / "environment.yml",
        }
        for module, (environment_reference, package_pin) in environment_checks.items():
            text = module.read_text(encoding="utf-8")
            self.assertIn(environment_reference, text)
            self.assertIn(CONTAINERS[module], text)
            self.assertIn('export TMPDIR="\\$PWD/.tmp"', text)
            self.assertIn('mkdir -p "\\$TMPDIR"', text)
            self.assertIn(package_pin, environment_paths[module].read_text(encoding="utf-8"))
            for forbidden in (
                "runtimeShellExports",
                "runtime.env_prefix",
                "${projectDir}",
                "/home/annan",
            ):
                self.assertNotIn(forbidden, text, f"{forbidden} in {module}")
        self.assertIn("base64 -d", FASTQC.read_text(encoding="utf-8"))
        self.assertNotIn("base64 --decode", FASTQC.read_text(encoding="utf-8"))

    def test_repository_qc_helper_and_multiqc_config_are_staged_inputs(self):
        samtools = SAMTOOLS_QC.read_text(encoding="utf-8")
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("path runtimeScripts, stageAs: 'tresflow/runtime/*'", samtools)
        self.assertIn('bash "${coreScriptsDir}/SamtoolsBamQc.sh"', samtools)
        self.assertIn("scripts/core_runtime/SamtoolsBamQc.sh", workflow)
        self.assertIn(
            "SAMTOOLS_BAM_QC(ch_bams_for_samtools_qc, samtoolsQcRuntimeScripts)",
            workflow,
        )
        self.assertIn('file("${projectDir}/assets/multiqc_config.yml")', workflow)
        self.assertIn("path(multiqc_config", MULTIQC.read_text(encoding="utf-8"))

    def test_samtools_wrapper_uses_path_and_preserves_commands(self):
        wrapper = (REPO / "scripts/core_runtime/SamtoolsBamQc.sh").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("SAMTOOLS_BIN", wrapper)
        for command in (
            "samtools quickcheck",
            "samtools flagstat",
            "samtools stats",
            "samtools idxstats",
        ):
            self.assertIn(command, wrapper)
        self.assertIn('idxstats_threads=$((threads > 0 ? threads - 1 : 0))', wrapper)

    def test_no_launch_time_host_tool_contract_remains(self):
        runtime = (REPO / "lib/RuntimeSupport.groovy").read_text(encoding="utf-8")
        exports = (REPO / "modules/local/runtime_support/main.nf").read_text(
            encoding="utf-8"
        )
        entry = (REPO / "main.nf").read_text(encoding="utf-8")
        self.assertNotIn("[name: 'python3', binary: 'python3']", runtime)
        self.assertNotIn("[name: 'samtools', binary: 'samtools']", runtime)
        self.assertNotIn("SAMTOOLS_BIN", runtime + exports)
        self.assertNotIn("PYTHON3_BIN", runtime + exports)
        self.assertNotIn("runtimeShellExports", exports)
        self.assertIn("writeCanonicalChromosomeContracts", runtime)
        self.assertNotIn("ProcessBuilder(command)", runtime)
        self.assertNotIn("runtimeToolPath", runtime)
        self.assertIn("runtimeSupport.writeCanonicalChromosomeContracts", entry)

        active_files = [
            WORKFLOW,
            REPO / "subworkflows/local/rna_core/main.nf",
            REPO / "subworkflows/local/dna_core/main.nf",
        ]
        active_graph = "\n".join(path.read_text(encoding="utf-8") for path in active_files)
        self.assertNotIn("modules/local/bam_coverage_dna", active_graph)
        self.assertNotIn("modules/local/mark_duplicates_dna", active_graph)

        active_tool_modules = [
            REPO / "modules/local/barcode_gate_metrics/main.nf",
            REPO / "modules/local/tres_report_html/main.nf",
            REPO / "modules/local/dual_tag_artifact_filter/main.nf",
            REPO / "modules/local/fq_to_sam/main.nf",
            REPO / "modules/local/split_dna_reads/main.nf",
            REPO / "modules/local/split_rna_reads/main.nf",
            REPO / "modules/local/tag_dna_cell_barcode/main.nf",
            REPO / "modules/local/tag_dna_modality/main.nf",
            REPO / "modules/local/tag_dna_sb/main.nf",
            REPO / "modules/local/tag_rna_cell_barcode/main.nf",
            REPO / "modules/local/tag_rna_sb/main.nf",
            REPO / "modules/local/tag_rna_umi/main.nf",
            REPO / "modules/local/trim_dna_fastqs/main.nf",
            REPO / "modules/local/trim_rna_fastqs/main.nf",
            REPO / "modules/local/rna_filtered_bam/main.nf",
            REPO / "modules/local/align_dna/main.nf",
            REPO / "modules/local/filter_canonical_dna_aligned_bam/main.nf",
            REPO / "modules/local/normalize_dna_markduplicates/main.nf",
            REPO / "modules/local/split_duplicates_dna/main.nf",
            REPO / "modules/nf-core/deeptools/bamcoverage/main.nf",
            REPO / "modules/nf-core/gatk4/markduplicates/main.nf",
            FASTQC,
            SAMTOOLS_QC,
            MULTIQC,
        ]
        for module in active_tool_modules:
            text = module.read_text(encoding="utf-8")
            self.assertIn("conda ", text, f"active task lacks Conda: {module}")
            self.assertNotIn(
                "runtimeShellExports", text, f"active task uses host exports: {module}"
            )
            self.assertNotIn("/home/annan", text, f"active task embeds a host path: {module}")

    def test_runtime_manifest_records_pins_hashes_licenses_and_platform(self):
        manifest = json.loads(
            (REPO / "modules/local/qc_processing/runtime-manifest.json").read_text()
        )
        self.assertEqual(manifest["supported_platforms"], ["linux/amd64"])
        packages = {item["name"]: item for item in manifest["conda_direct_packages"]}
        self.assertEqual(set(packages), {"fastqc", "samtools", "multiqc"})
        self.assertEqual(packages["fastqc"]["version"], "0.12.1")
        self.assertEqual(packages["samtools"]["version"], "1.23.1")
        self.assertEqual(packages["multiqc"]["version"], "1.33")
        for package in packages.values():
            self.assertRegex(package["sha256"], r"^[0-9a-f]{64}$")
            self.assertTrue(package["license"])
        for container in manifest["containers"].values():
            self.assertIn("@sha256:", container["reference"])
            self.assertEqual(container["platform"], "linux/amd64")

    def test_apptainer_profile_and_real_isolated_harness_are_explicit(self):
        config = (REPO / "tests/phase6_qc.config").read_text(encoding="utf-8")
        harness = (REPO / "tests/phase6_qc.nf").read_text(encoding="utf-8")
        self.assertIn("phase6_qc_apptainer_static", config)
        self.assertIn("apptainer.enabled = true", config)
        self.assertIn("apptainer.autoMounts = true", config)
        self.assertNotIn("runtime_env_prefix", harness)
        self.assertNotIn("runtime_tmpdir", harness)
        for process in ("FASTQC", "SAMTOOLS_BAM_QC", "TRES_REPORT_HTML", "MULTIQC"):
            self.assertIn(f"{process}(", harness)
        self.assertIn("I2.fastq", harness)
        self.assertIn("dual_tag_artifact_filter.cutadapt.json", harness)
        self.assertIn("DuplicateMetrics.txt", harness)
        self.assertIn("Solo.outGeneFull/Summary.csv", harness)


if __name__ == "__main__":
    unittest.main()
