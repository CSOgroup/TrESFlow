from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def read(relative):
    return (REPO / relative).read_text(encoding="utf-8")


def test_supported_samplesheets_do_not_declare_a_runtime_environment():
    for relative in (
        "assets/samplesheet.example.yaml",
        "assets/samplesheet.real.example.yaml",
        "assets/samplesheet.template.yaml",
    ):
        text = read(relative)
        assert "\nruntime:" not in text
        assert "env_prefix:" not in text
        assert "tmpdir:" not in text


def test_legacy_samplesheet_runtime_fields_are_compatibility_only():
    parser = read("lib/SamplesheetParser.groovy")
    entry = read("main.nf")

    assert "resolveDeprecatedRuntimeFields" in parser
    assert "deprecated_runtime_fields" in parser
    assert "Deprecated samplesheet runtime fields are ignored" in entry
    for forbidden in (
        "row.runtime_env_prefix",
        "row.runtime_tmpdir",
        "runtimeConfig",
        "runtimeParams",
        "validateRuntimeContract",
    ):
        assert forbidden not in parser + entry


def test_host_preflight_exports_and_dead_modules_are_removed():
    runtime = read("lib/RuntimeSupport.groovy")
    task_support = read("modules/local/runtime_support/main.nf")
    workflow = read("workflows/treseq.nf")

    for token in (
        "runtimeShellExports",
        "runtimeCoreScriptsDir",
        "runtimeToolPath",
        "runtimeEnvPrefix",
        "runtimeTmpdir",
        "PYTHON3_BIN",
        "SAMTOOLS_BIN",
        "CODON_BIN",
        "CODON_HOME",
        "RUNTIME_ENV_PREFIX",
        "RUNTIME_BIN_DIR",
    ):
        assert token not in runtime + task_support + workflow
    assert "runtimeOutdir" in task_support
    assert not (REPO / "bin/check_codon_seq_host.sh").exists()
    assert not (REPO / "bin/resolve_canonical_chromosomes.py").exists()
    assert not (REPO / "modules/local/bam_coverage_dna/main.nf").exists()
    assert not (REPO / "modules/local/mark_duplicates_dna/main.nf").exists()


def test_active_canonical_filter_uses_samtools_from_path():
    helper = read("scripts/core_runtime/FilterCanonicalBam.sh")
    rna = read("scripts/core_runtime/RNA_FILTERED_BAM.sh")
    dna = read("modules/local/filter_canonical_dna_aligned_bam/main.nf")
    markeddup = read("modules/local/normalize_dna_markduplicates/main.nf")

    assert "SAMTOOLS_BIN" not in helper + rna + dna + markeddup
    assert "samtools view" in helper
    assert 'bash "${script_dir}/FilterCanonicalBam.sh"' in rna


def test_documentation_describes_process_environments_and_work_dir():
    docs = "\n".join(
        read(relative)
        for relative in (
            "README.md",
            "docs/usage.md",
            "docs/output.md",
            "docs/architecture/implemented_pipeline.md",
        )
    )
    assert "-with-conda" in docs
    assert "-work-dir" in docs
    assert "accepted but ignored" in docs
    assert "Create and activate a conda/mamba/micromamba environment" not in docs
