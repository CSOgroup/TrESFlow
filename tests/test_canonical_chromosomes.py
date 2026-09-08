import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNTIME_SUPPORT = REPO_ROOT / "lib/RuntimeSupport.groovy"
FIXTURES = REPO_ROOT / "tests/fixtures/canonical_chromosomes"

GROOVY_RESOLVE = r'''
def runtime = new GroovyClassLoader().parseClass(new File(args[0]))
def command = args[1]
try {
    if (command == 'chrom-sizes') {
        def resolved = runtime.resolveCanonicalEntries(
            runtime.readChromSizes(new File(args[2])),
            args[2]
        )
        println groovy.json.JsonOutput.toJson(resolved)
    }
    else if (command == 'write') {
        def outdir = args[2]
        def rnaSizes = args[3]
        def dnaPrefix = args[4]
        def dnaSizes = args[5]
        def references = [
            rna_chrom_sizes  : rnaSizes == '-' ? null : rnaSizes,
            dna_bwa_reference: dnaPrefix == '-' ? null : dnaPrefix,
            dna_chrom_sizes  : dnaSizes == '-' ? null : dnaSizes,
        ]
        def modalities = [rna: rnaSizes != '-', dna: dnaPrefix != '-']
        println groovy.json.JsonOutput.toJson(
            runtime.writeCanonicalChromosomeContracts(outdir, references, modalities)
        )
    }
    else {
        throw new IllegalArgumentException("Unknown test command: ${command}")
    }
}
catch (Exception error) {
    System.err.println(error.message)
    System.exit(2)
}
'''


def find_nextflow_jar():
    candidates = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.extend(Path(conda_prefix).glob("share/nextflow/dist/*/nextflow-*-one.jar"))
    candidates.extend(
        Path("/home/annan/micromamba/envs/tres/share/nextflow/dist").glob(
            "*/nextflow-*-one.jar"
        )
    )
    candidates.extend(Path.home().glob(".nextflow/framework/*/nextflow-*-one.jar"))
    return sorted(candidates)[-1] if candidates else None


def run_runtime(*args, check=True):
    java = shutil.which("java")
    jar = find_nextflow_jar()
    if not java or jar is None:
        pytest.skip("Java and a local Nextflow distribution are required")
    return subprocess.run(
        [
            java,
            "-cp",
            str(jar),
            "groovy.ui.GroovyMain",
            "-e",
            GROOVY_RESOLVE,
            "--",
            str(RUNTIME_SUPPORT),
            *map(str, args),
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=check,
    )


@pytest.mark.parametrize(
    ("fixture", "expected_style", "expected_names"),
    [
        (
            "hg19_ucsc.chrom.sizes",
            "ucsc",
            ["chr1", "chrX", "chrY", "chrM"],
        ),
        (
            "hg38_ucsc.chrom.sizes",
            "ucsc",
            ["chr1", "chrX", "chrY", "chrM"],
        ),
        (
            "hg38_ensembl.chrom.sizes",
            "ensembl",
            ["1", "X", "Y", "MT"],
        ),
    ],
)
def test_groovy_resolver_keeps_exact_names_lengths_and_reference_order(
    fixture, expected_style, expected_names
):
    result = run_runtime("chrom-sizes", FIXTURES / fixture)
    resolved = json.loads(result.stdout)

    assert resolved["style"] == expected_style
    assert [entry["name"] for entry in resolved["entries"]] == expected_names
    source_lengths = {
        line.split()[0]: int(line.split()[1])
        for line in (FIXTURES / fixture).read_text().splitlines()
    }
    assert {
        entry["name"]: entry["length"] for entry in resolved["entries"]
    } == {name: source_lengths[name] for name in expected_names}


def test_bwa_ann_and_chrom_sizes_write_matching_exact_contracts(tmp_path):
    result = run_runtime(
        "write",
        tmp_path,
        FIXTURES / "hg19_ucsc.chrom.sizes",
        FIXTURES / "hg38_ucsc.fa",
        FIXTURES / "hg38_ucsc.chrom.sizes",
    )
    contracts = json.loads(result.stdout)

    assert contracts["rna"]["style"] == "ucsc"
    assert contracts["dna"]["style"] == "ucsc"
    assert contracts["rna"]["contigs"] == ["chr1", "chrX", "chrY", "chrM"]
    assert contracts["dna"]["contigs"] == ["chr1", "chrX", "chrY", "chrM"]
    derived = tmp_path / "pipeline_info/derived_contract"
    assert (derived / "rna_canonical_chromosomes.txt").read_bytes() == (
        b"chr1\nchrX\nchrY\nchrM\n"
    )
    assert (derived / "dna_canonical_chromosomes.chrom.sizes").read_bytes() == (
        b"chr1\t248956422\nchrX\t156040895\nchrY\t57227415\nchrM\t16569\n"
    )


def test_bwa_and_explicit_chromosome_sizes_mismatch_fails(tmp_path):
    mismatched = tmp_path / "mismatched.chrom.sizes"
    mismatched.write_text(
        "chr1\t248956421\nchrX\t156040895\nchrY\t57227415\nchrM\t16569\n"
    )
    result = run_runtime(
        "write",
        tmp_path / "output",
        "-",
        FIXTURES / "hg38_ucsc.fa",
        mismatched,
        check=False,
    )

    assert result.returncode == 2
    assert "Canonical chromosome names or lengths disagree" in result.stderr


def test_mixed_conventions_and_missing_anchors_fail_clearly(tmp_path):
    mixed = tmp_path / "mixed.chrom.sizes"
    mixed.write_text("chr1\t100\nX\t90\nY\t80\nMT\t70\n")
    result = run_runtime("chrom-sizes", mixed, check=False)
    assert result.returncode == 2
    assert "both UCSC-style" in result.stderr
    assert "Ensembl-style" in result.stderr

    missing = tmp_path / "missing.chrom.sizes"
    missing.write_text("chr1\t100\nchrX\t90\nchrY\t80\n")
    result = run_runtime("chrom-sizes", missing, check=False)
    assert result.returncode == 2
    assert "missing: chrM" in result.stderr


def test_resolver_is_pure_groovy_and_runtime_environment_plumbing_is_gone():
    runtime = RUNTIME_SUPPORT.read_text()
    entry = (REPO_ROOT / "main.nf").read_text()
    task_support = (REPO_ROOT / "modules/local/runtime_support/main.nf").read_text()

    for forbidden in (
        "resolve_canonical_chromosomes.py",
        "runtimeToolPath",
        "validateRuntimeContract",
        "runtimeEnvPrefix",
        "runtimeTmpdir",
        "PYTHON3_BIN",
        "RUNTIME_ENV_PREFIX",
        "runtimeShellExports",
    ):
        assert forbidden not in runtime + entry + task_support
    assert "readChromSizes" in runtime
    assert "readBwaAnn" in runtime
    assert "runtimeOutdir" in task_support


def test_rna_coverage_does_not_force_ucsc_prefix():
    coverage_script = (
        REPO_ROOT / "scripts/core_runtime/RNA_COVERAGE.sh"
    ).read_text(encoding="utf-8")
    workflow_text = (
        REPO_ROOT / "subworkflows/local/rna_core/main.nf"
    ).read_text(encoding="utf-8")

    assert "--outWigReferencesPrefix" not in coverage_script
    assert "meta.canonical_chrom_sizes" in workflow_text
