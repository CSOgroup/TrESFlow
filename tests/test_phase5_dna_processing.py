import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[1]
ALIGN = REPO / "modules/local/align_dna/main.nf"
FILTER = REPO / "modules/local/filter_canonical_dna_aligned_bam/main.nf"
NORMALIZE = REPO / "modules/local/normalize_dna_markduplicates/main.nf"
SPLIT = REPO / "modules/local/split_duplicates_dna/main.nf"
GATK = REPO / "modules/nf-core/gatk4/markduplicates/main.nf"
COVERAGE = REPO / "modules/nf-core/deeptools/bamcoverage/main.nf"
MODULES = (ALIGN, FILTER, NORMALIZE, SPLIT, GATK, COVERAGE)
CONTAINERS = {
    ALIGN: (
        "community.wave.seqera.io/library/bwa-mem2_samtools@"
        "sha256:ce8cbf5cc21c690c8c2994d9bbb409b9313c47f38c5452a5fe7dec3402eff9c8"
    ),
    FILTER: (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    NORMALIZE: (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    SPLIT: (
        "community.wave.seqera.io/library/samtools@"
        "sha256:2ee310db4ac650bc54c16dc9d28151d973e2ffed0ca878de8fc8e70e820ffe34"
    ),
    GATK: (
        "community.wave.seqera.io/library/gatk4_samtools@"
        "sha256:ca703ab322c8cf829f3987275648ab98fec87b4a340b5cbe154a49cb9a4f41a0"
    ),
    COVERAGE: (
        "community.wave.seqera.io/library/deeptools_samtools@"
        "sha256:9c482aa632f9a30dcd09c38bed01ac89c54b345fd1ec666deb3ec9c595b258c9"
    ),
}


class Phase5DnaProcessingArchitectureTests(unittest.TestCase):
    def test_active_processes_use_exact_environments_and_containers(self):
        environments = {
            ALIGN: "environment-align.yml",
            FILTER: "environment-samtools.yml",
            NORMALIZE: "environment-samtools.yml",
            SPLIT: "environment-samtools.yml",
            GATK: "environment.yml",
            COVERAGE: "environment.yml",
        }
        for module in MODULES:
            text = module.read_text(encoding="utf-8")
            self.assertIn(environments[module], text)
            self.assertIn(CONTAINERS[module], text)
            self.assertIn('export TMPDIR="\\$PWD/.tmp"', text)
            self.assertIn('mkdir -p "\\$TMPDIR"', text)
            for forbidden in (
                "runtimeShellExports",
                "runtime.env_prefix",
                "BWA_MEM2_BIN",
                "GATK_BIN",
                "BAMCOVERAGE_BIN",
                "${projectDir}",
                "/home/annan",
            ):
                self.assertNotIn(forbidden, text, f"{forbidden} in {module}")

    def test_repository_helpers_and_external_bwa_index_are_explicit_inputs(self):
        align = ALIGN.read_text(encoding="utf-8")
        self.assertIn("path(bwaIndexFiles, stageAs: 'bwa_index/*')", align)
        self.assertIn("path runtimeScripts, stageAs: 'tresflow/runtime/*'", align)
        for module in (FILTER, NORMALIZE, SPLIT):
            self.assertIn(
                "path runtimeScripts, stageAs: 'tresflow/runtime/*'",
                module.read_text(encoding="utf-8"),
            )
        self.assertIn(
            'chmod a+r "${splitName}.bam" "${splitName}.bam.bai"',
            FILTER.read_text(encoding="utf-8"),
        )
        self.assertIn(
            'chmod a+r "${splitName}_MarkedDup.bam" "${splitName}_MarkedDup.bam.bai"',
            NORMALIZE.read_text(encoding="utf-8"),
        )

        workflow = (REPO / "subworkflows/local/dna_core/main.nf").read_text(
            encoding="utf-8"
        )
        for helper in ("AlignDNA.sh", "FilterCanonicalBam.sh", "SplitDuplicatesDNA.sh"):
            self.assertIn(helper, workflow)
        for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
            self.assertIn(suffix, workflow)
        self.assertIn("file(\"${bwaReference}${suffix}\", checkIfExists: true)", workflow)
        self.assertNotIn("genomeGenerate", workflow)

    def test_wrappers_invoke_path_tools_and_preserve_alignment_operations(self):
        align = (REPO / "scripts/core_runtime/AlignDNA.sh").read_text(encoding="utf-8")
        split = (REPO / "scripts/core_runtime/SplitDuplicatesDNA.sh").read_text(
            encoding="utf-8"
        )
        for forbidden in ("BWA_MEM2_BIN", "SAMTOOLS_BIN"):
            self.assertNotIn(forbidden, align)
            self.assertNotIn(forbidden, split)
        self.assertIn("bwa-mem2 mem -t ${threads} -C -o", align)
        self.assertIn("--require-flags 0x2", align)
        self.assertIn("-L ${blacklist_bed}", align)
        self.assertIn("samtools sort --threads ${sort_threads} -m ${sort_mem} -l 0 -n", align)
        self.assertIn("samtools sort -@ ${sort_threads} -m ${sort_mem}", align)
        self.assertIn("--exclude-flags 0x400", split)
        self.assertIn("samtools index", split)
        self.assertIn("samtools view --threads \"${threads}\" -c -F 4", split)

    def test_markduplicates_arguments_and_output_contract_are_unchanged(self):
        config = (REPO / "conf/modules.config").read_text(encoding="utf-8")
        for argument in (
            "--REMOVE_DUPLICATES false",
            "--BARCODE_TAG CB",
            "--READ_NAME_REGEX",
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.aviti_optical_duplicate_distance}",
            "--CREATE_INDEX true",
            "--MAX_RECORDS_IN_RAM 10000000",
        ):
            self.assertIn(argument, config)
        self.assertIn(
            "aviti_optical_duplicate_distance = 10",
            (REPO / "nextflow.config").read_text(encoding="utf-8"),
        )
        gatk = GATK.read_text(encoding="utf-8")
        self.assertIn("gatk --java-options", gatk)
        self.assertIn("MarkDuplicates", gatk)
        self.assertIn("--TMP_DIR .", gatk)

    def test_nodup_bams_are_the_only_runtime_coverage_source(self):
        workflow = (REPO / "subworkflows/local/dna_core/main.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("ch_nodup_for_coverage = SPLIT_DUPLICATES_DNA.out.bam", workflow)
        self.assertIn(".join(SPLIT_DUPLICATES_DNA.out.bai)", workflow)
        self.assertIn(".join(SPLIT_DUPLICATES_DNA.out.mapped_reads)", workflow)
        self.assertIn("hasMappedDnaReads(mappedReadsFile)", workflow)
        self.assertIn("DEEPTOOLS_BAMCOVERAGE(\n        ch_nodup_for_coverage", workflow)
        coverage = COVERAGE.read_text(encoding="utf-8")
        self.assertIn("bamCoverage \\", coverage)
        self.assertIn("--bam $input_out", coverage)

    def test_exact_conda_packages_and_runtime_manifest(self):
        runtime_root = REPO / "modules/local/dna_processing"
        align_env = (runtime_root / "environment-align.yml").read_text(encoding="utf-8")
        samtools_env = (runtime_root / "environment-samtools.yml").read_text(
            encoding="utf-8"
        )
        gatk_env = (GATK.parent / "environment.yml").read_text(encoding="utf-8")
        coverage_env = (COVERAGE.parent / "environment.yml").read_text(encoding="utf-8")
        self.assertIn("bioconda::bwa-mem2=2.2.1=he70b90d_8", align_env)
        self.assertIn("bioconda::samtools=1.23.1=ha83d96e_0", align_env)
        self.assertIn("bioconda::samtools=1.23.1=ha83d96e_0", samtools_env)
        self.assertIn("bioconda::gatk4=4.6.2.0=py310hdfd78af_0", gatk_env)
        self.assertIn("bioconda::samtools=1.23.1=ha83d96e_0", gatk_env)
        self.assertIn("bioconda::deeptools=3.5.5=pyhdfd78af_0", coverage_env)
        self.assertIn("bioconda::samtools=1.23.1=ha83d96e_0", coverage_env)

        manifest = json.loads((runtime_root / "runtime-manifest.json").read_text())
        self.assertEqual(manifest["supported_platforms"], ["linux/amd64"])
        self.assertEqual(
            {item["name"] for item in manifest["conda_direct_packages"]},
            {"bwa-mem2", "samtools", "gatk4", "deeptools"},
        )
        self.assertEqual(
            {
                item["name"]: item["license"]
                for item in manifest["conda_direct_packages"]
            },
            {
                "bwa-mem2": "MIT",
                "samtools": "MIT/Expat",
                "gatk4": "BSD-3-Clause",
                "deeptools": "MIT (deeptools/cm.py: BSD)",
            },
        )
        for container in manifest["containers"].values():
            self.assertIn("@sha256:", container["reference"])
            self.assertEqual(container["platform"], "linux/amd64")
        for package in manifest["conda_direct_packages"]:
            self.assertRegex(package["sha256"], r"^[0-9a-f]{64}$")

    def test_global_host_contract_retires_only_completed_family_tools(self):
        runtime = (REPO / "lib/RuntimeSupport.groovy").read_text(encoding="utf-8")
        exports = (REPO / "modules/local/runtime_support/main.nf").read_text(
            encoding="utf-8"
        )
        for retired in (
            "[name: 'bwa-mem2', binary: 'bwa-mem2']",
            "[name: 'bamCoverage', binary: 'bamCoverage']",
            "[name: 'gatk', binary: 'gatk']",
            "BWA_MEM2_BIN",
            "BAMCOVERAGE_BIN",
            "GATK_BIN",
        ):
            self.assertNotIn(retired, runtime + exports)
        for retained in ("python3", "samtools", "PYTHON3_BIN", "SAMTOOLS_BIN"):
            self.assertIn(retained, runtime + exports)

    def test_obsolete_dna_modules_are_not_in_the_active_graph(self):
        workflow = (REPO / "subworkflows/local/dna_core/main.nf").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("modules/local/mark_duplicates_dna", workflow)
        self.assertNotIn("modules/local/bam_coverage_dna", workflow)
        self.assertIn("modules/nf-core/gatk4/markduplicates", workflow)
        self.assertIn("modules/nf-core/deeptools/bamcoverage", workflow)

    def test_apptainer_profile_and_isolated_harness_are_explicit(self):
        config = (REPO / "tests/phase5_dna_processing.config").read_text(
            encoding="utf-8"
        )
        harness = (REPO / "tests/phase5_dna_processing.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("phase5_dna_apptainer_static", config)
        self.assertIn("apptainer.enabled = true", config)
        self.assertIn("apptainer.autoMounts = true", config)
        self.assertIn("phase5-must-not-be-used", harness)
        self.assertIn("empty_markeddup_bam", harness)
        self.assertIn("hasMappedDnaReads(mappedReads)", harness)
        for process_name in (
            "ALIGN_DNA",
            "FILTER_CANONICAL_DNA_ALIGNED_BAM",
            "GATK4_MARKDUPLICATES",
            "NORMALIZE_DNA_MARKDUPLICATES",
            "SPLIT_DUPLICATES_DNA",
            "DEEPTOOLS_BAMCOVERAGE",
        ):
            self.assertIn(f"{process_name}(", harness)

    def test_semantic_validator_covers_phase0_and_runtime_sources(self):
        validator = (REPO / "tests/validate_phase5_dna_processing.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("normalizer.capture_contract", validator)
        self.assertIn("assert_bam_invariants", validator)
        self.assertIn("assert_duplicate_metrics", validator)
        self.assertIn("_NoDup.bam", validator)
        self.assertIn("phase5_empty", validator)
        self.assertIn("docker run ", validator)
        self.assertIn("micromamba activate ", validator)


if __name__ == "__main__":
    unittest.main()
