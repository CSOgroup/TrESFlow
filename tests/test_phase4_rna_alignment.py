import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[1]
STARSOLO_MODULE = REPO / "modules/local/rna_starsolo_align/main.nf"
FILTER_MODULE = REPO / "modules/local/rna_filtered_bam/main.nf"
COVERAGE_MODULE = REPO / "modules/local/rna_coverage/main.nf"
MODULES = (STARSOLO_MODULE, FILTER_MODULE, COVERAGE_MODULE)
STAR_CONTAINER = (
    "quay.io/biocontainers/star@"
    "sha256:d5139d9e29bc8a871a020ca5b21bc0599fac428bd061d4b13d863f4cb4400d31"
)
FILTER_CONTAINER = (
    "community.wave.seqera.io/library/samtools_python@"
    "sha256:74535d380b6c327aa8a82ad941f00d900d26f1a74217e82679f9d64b1b9e28d3"
)
COVERAGE_CONTAINER = (
    "community.wave.seqera.io/library/star_ucsc-bedgraphtobigwig@"
    "sha256:133ac55ecc30285f3e8b0efa7e3577efd6184354641098d0a290c66268038d88"
)


class Phase4RnaAlignmentArchitectureTests(unittest.TestCase):
    def test_active_processes_use_exact_family_environments_and_containers(self):
        expected = {
            STARSOLO_MODULE: ("environment-star.yml", STAR_CONTAINER),
            FILTER_MODULE: ("environment-filter.yml", FILTER_CONTAINER),
            COVERAGE_MODULE: ("environment-coverage.yml", COVERAGE_CONTAINER),
        }
        for module, (environment, container) in expected.items():
            text = module.read_text(encoding="utf-8")
            self.assertIn(environment, text)
            self.assertIn(container, text)
            self.assertIn("label 'rna_alignment'", text)
            self.assertIn("path runtimeScripts, stageAs: 'tresflow/runtime/*'", text)
            self.assertIn('export TMPDIR="\\$PWD/.tmp"', text)
            self.assertIn('mkdir -p "\\$TMPDIR"', text)
            for forbidden in (
                "runtimeShellExports",
                "runtimeCoreScriptsDir",
                "runtime.env_prefix",
                "${projectDir}",
                "/home/annan",
            ):
                self.assertNotIn(forbidden, text, f"{forbidden} in {module}")

    def test_external_references_are_declared_as_task_path_inputs(self):
        starsolo = STARSOLO_MODULE.read_text(encoding="utf-8")
        coverage = COVERAGE_MODULE.read_text(encoding="utf-8")
        self.assertIn("path(starIndexDir)", starsolo)
        self.assertIn("path(starIndexDir)", coverage)
        self.assertIn("path(chromSizes)", coverage)
        self.assertIn('chmod -R a+rX "${splitName}.Solo.outGeneFull"', starsolo)
        for text in (starsolo, coverage):
            self.assertNotIn("COPY", text)
            self.assertNotIn("genomeGenerate", text)

    def test_repository_scripts_are_staged_by_the_production_subworkflow(self):
        workflow = (REPO / "subworkflows/local/rna_core/main.nf").read_text(
            encoding="utf-8"
        )
        for name in (
            "RNA_STARSOLO_ALIGN.sh",
            "RNA_FILTERED_BAM.sh",
            "SummarizeRnaRetention.py",
            "FilterCanonicalBam.sh",
            "RNA_COVERAGE.sh",
        ):
            self.assertIn(name, workflow)
        self.assertIn(
            "RNA_STARSOLO_ALIGN(ch_starsolo_input, starsoloRuntimeScripts)", workflow
        )
        self.assertIn(
            "RNA_FILTERED_BAM(ch_filtered_bam_input, filteredBamRuntimeScripts)",
            workflow,
        )
        self.assertIn("RNA_COVERAGE(ch_coverage_input, coverageRuntimeScripts)", workflow)
        self.assertIn("file(meta.rna_star_index_dir)", workflow)
        self.assertIn("file(meta.canonical_chrom_sizes)", workflow)

    def test_rna_runtime_wrappers_use_normal_path_without_host_overrides(self):
        starsolo = (REPO / "scripts/core_runtime/RNA_STARSOLO_ALIGN.sh").read_text()
        filtering = (REPO / "scripts/core_runtime/RNA_FILTERED_BAM.sh").read_text()
        coverage = (REPO / "scripts/core_runtime/RNA_COVERAGE.sh").read_text()

        self.assertIn("\nSTAR \\\n", starsolo)
        self.assertNotIn("STAR_BIN", starsolo)
        self.assertIn('samtools view "${INBAM}"', filtering)
        self.assertIn('python3 "${script_dir}/SummarizeRnaRetention.py"', filtering)
        self.assertIn("env -u SAMTOOLS_BIN", filtering)
        self.assertIn(
            'chmod a+r "${splitName}.filtered_cells.bam"',
            FILTER_MODULE.read_text(encoding="utf-8"),
        )
        self.assertNotIn("${SAMTOOLS_BIN}", filtering)
        self.assertNotIn("PYTHON3_BIN", filtering)
        self.assertIn("\n    STAR \\\n", coverage)
        self.assertIn("\n    bedGraphToBigWig \\\n", coverage)
        self.assertNotIn("STAR_BIN", coverage)
        self.assertNotIn("BEDGRAPH_TO_BIGWIG_BIN", coverage)

    def test_star_arguments_and_coverage_transform_remain_present(self):
        starsolo = (REPO / "scripts/core_runtime/RNA_STARSOLO_ALIGN.sh").read_text()
        coverage = (REPO / "scripts/core_runtime/RNA_COVERAGE.sh").read_text()
        for argument in (
            "--twopassMode Basic",
            "--soloFeatures GeneFull",
            "--soloStrand Forward",
            "--soloMultiMappers EM",
            "--soloUMIdedup 1MM_CR",
            "--soloUMIfiltering MultiGeneUMI_CR",
            "--outSAMtype BAM SortedByCoordinate",
            "--outBAMcompression 0",
        ):
            self.assertIn(argument, starsolo)
        for operation in (
            "--runMode inputAlignmentsFromBAM",
            "--outWigStrand Stranded",
            "--outWigStrand Unstranded",
            "--outWigNorm RPM",
            "sort -k1,1 -k2,2n",
        ):
            self.assertIn(operation, coverage)

    def test_exact_conda_packages_and_runtime_manifest(self):
        runtime_root = REPO / "modules/local/rna_alignment"
        expected = {
            "environment-star.yml": ("bioconda::star=2.7.11b=h5ca1c30_7",),
            "environment-filter.yml": (
                "conda-forge::python=3.12.13=hd63d673_0_cpython",
                "bioconda::samtools=1.23.1=ha83d96e_0",
            ),
            "environment-coverage.yml": (
                "bioconda::star=2.7.11b=h5ca1c30_7",
                "bioconda::ucsc-bedgraphtobigwig=482=hdc0a859_0",
            ),
        }
        for filename, requirements in expected.items():
            environment = (runtime_root / filename).read_text(encoding="utf-8")
            for requirement in requirements:
                self.assertIn(requirement, environment)

        manifest = json.loads((runtime_root / "runtime-manifest.json").read_text())
        self.assertEqual(manifest["supported_platforms"], ["linux/amd64"])
        self.assertEqual(
            manifest["containers"]["starsolo_alignment"]["reference"],
            STAR_CONTAINER,
        )
        self.assertEqual(
            manifest["containers"]["rna_filtering"]["reference"], FILTER_CONTAINER
        )
        self.assertEqual(
            manifest["containers"]["rna_coverage"]["reference"],
            COVERAGE_CONTAINER,
        )
        self.assertEqual(
            {package["name"] for package in manifest["conda_direct_packages"]},
            {"python", "samtools", "star", "ucsc-bedgraphtobigwig"},
        )
        for package in manifest["conda_direct_packages"]:
            self.assertRegex(package["sha256"], r"^[0-9a-f]{64}$")

    def test_global_host_contract_retires_only_migrated_tools(self):
        runtime = (REPO / "lib/RuntimeSupport.groovy").read_text(encoding="utf-8")
        exports = (REPO / "modules/local/runtime_support/main.nf").read_text(
            encoding="utf-8"
        )
        for retired in (
            "[name: 'STAR', binary: 'STAR']",
            "[name: 'bedGraphToBigWig', binary: 'bedGraphToBigWig']",
            "STAR_BIN",
            "BEDGRAPH_TO_BIGWIG_BIN",
        ):
            self.assertNotIn(retired, runtime + exports)
        for retained in ("python3", "samtools"):
            self.assertIn(retained, runtime)

    def test_apptainer_profile_reuses_digest_pinned_oci_images(self):
        config = (REPO / "tests/phase4_rna_alignment.config").read_text(
            encoding="utf-8"
        )
        self.assertIn("phase4_rna_apptainer_static", config)
        self.assertIn("apptainer.enabled = true", config)
        self.assertIn("apptainer.autoMounts = true", config)
        for module in MODULES:
            self.assertIn("@sha256:", module.read_text(encoding="utf-8"))

    def test_isolated_real_harness_poisoned_the_legacy_host_prefix(self):
        harness = (REPO / "tests/phase4_rna_alignment.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("phase4-must-not-be-used", harness)
        for process_name in (
            "RNA_STARSOLO_ALIGN",
            "RNA_FILTERED_BAM",
            "RNA_COVERAGE",
        ):
            self.assertIn(f"{process_name}(", harness)
        self.assertIn("file(params.star_index, checkIfExists: true)", harness)

    def test_isolated_validator_uses_semantic_phase0_decoders(self):
        validator = (REPO / "tests/validate_phase4_rna_alignment.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("normalizer.capture_contract", validator)
        self.assertIn('len(actual["bams"]) != 1', validator)
        self.assertIn('len(actual["bigwigs"]) != 3', validator)
        self.assertIn('len(actual["star_matrices"]) != 3', validator)


if __name__ == "__main__":
    unittest.main()
