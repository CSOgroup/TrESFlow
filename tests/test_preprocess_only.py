"""Mode-aware parser and real Nextflow orchestration tests.

The test profile mocks biological stages; FastQC, MultiQC, compression, metric
sidecars, publication and cleanup execute normally. The separate real-toolchain
smoke test runs every preprocessing stage with Codon/Seq and Trim Galore.
"""
import copy
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import shlex
import subprocess
from pathlib import Path

import pytest
import yaml

from test_multi_fastq_inputs import GROOVY_PARSE, PARSER, find_nextflow_jar

REPO = Path(__file__).resolve().parents[1]
RUNTIME = Path(os.environ.get("TRESFLOW_TEST_RUNTIME", "/home/annan/micromamba/envs/tres"))
FIXTURES = {
    "rna": REPO / "tests/samplesheets/multi_rna_comma.yaml",
    "single": REPO / "tests/samplesheets/multi_dna_single_sequence.yaml",
    "dual": REPO / "tests/samplesheets/multi_dna_dual_sequence.yaml",
    "combined": REPO / "tests/samplesheets/group_specific_modalities.yaml",
}
ALIGNMENT_TOOLS = {"STAR", "bwa-mem2", "samtools", "gatk", "bamCoverage", "bedGraphToBigWig"}
DISABLED = {"FQ_TO_SAM", "RNA_STARSOLO_ALIGN", "RNA_FILTERED_BAM", "RNA_COVERAGE", "ALIGN_DNA", "FILTER_CANONICAL_DNA_ALIGNED_BAM", "GATK4_MARKDUPLICATES", "NORMALIZE_DNA_MARKDUPLICATES", "SPLIT_DUPLICATES_DNA", "DEEPTOOLS_BAMCOVERAGE", "SAMTOOLS_BAM_QC"}


def fixture(kind):
    path = FIXTURES[kind]
    data = yaml.safe_load(path.read_text())
    for key, value in data["references"].items():
        if key not in {"species", "dna_effective_genome_size"}:
            data["references"][key] = str((path.parent / value).resolve())
    for sample in data["samples"].values():
        for modality in ("rna", "dna"):
            if modality not in sample:
                continue
            for role, paths in sample[modality]["reads"].items():
                paths = paths if isinstance(paths, list) else [p.strip() for p in paths.split(",")]
                sample[modality]["reads"][role] = [str((path.parent / p).resolve()) for p in paths]
    data["runtime"] = {"env_prefix": str(RUNTIME)}
    return data


def extend_rna_fixture(root, data):
    # Keep one long pair per technical chunk; shorter pairs exercise export
    # loss accounting. Barcode/UMI positions and index reads stay identical.
    for sample in data["samples"].values():
        if "rna" not in sample:
            continue
        for role in ("r1", "r2"):
            paths = []
            for index, path in enumerate(sample["rna"]["reads"][role]):
                lines = Path(path).read_text().splitlines()
                lines[1] += "GCTACGTCAGTCCGATGACCTGAGTCGACTGA" * 2
                lines[3] += "F" * 64
                destination = root / f"fixture_rna_{role}_{index}.fastq"
                destination.write_text("\n".join(lines) + "\n")
                paths.append(str(destination))
            sample["rna"]["reads"][role] = paths
    return data


def assert_rna_conversion(exported, original):
    for name, raw in exported.items():
        if not name.startswith("rna_") or "_R1." not in name:
            continue
        before1, before2 = original[name].decode().splitlines(), original[name.replace("_R1.", "_R2.")].decode().splitlines()
        one, two = raw.decode().splitlines(), exported[name.replace("_R1.", "_R2.")].decode().splitlines()
        kept = [offset for offset in range(0, len(before1), 4) if len(before1[offset + 1]) >= 20 and len(before2[offset + 1]) >= 34]
        assert len(one) == len(two) == len(kept) * 4
        for after, before in enumerate(kept):
            offset = after * 4
            tags = dict(token.split(":Z:", 1) for token in before1[before].split()[1:] if ":Z:" in token)
            pool = "P" + hashlib.sha256(json.dumps([tags["CB"], tags["SB"]], separators=(",", ":")).encode()).hexdigest()
            primary = before1[before].split()[0] + f"_{pool}_{tags['UM']}"
            assert one[offset].split()[0] == two[offset].split()[0] == primary
            assert one[offset + 1:offset + 4] == before1[before + 1:before + 4]
            assert two[offset + 1:offset + 4] == [before2[before + 1][14:], "+", before2[before + 3][14:]]


def reference_free(data):
    data = copy.deepcopy(data)
    data["references"] = {"ligation_barcode_whitelist": data["references"]["ligation_barcode_whitelist"]}
    return data


def write_sheet(root, data):
    root.mkdir(parents=True, exist_ok=True)
    path = root / "samples.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def parse_contract(root, data, preprocess=True):
    jar = find_nextflow_jar()
    if jar is None or not shutil.which("java"):
        pytest.skip("Java/local Nextflow jar unavailable for parser checks")
    script = GROOVY_PARSE.replace("barcode_defaults: defaults]", "barcode_defaults: defaults, preprocess_only: args[3].toBoolean()]")
    # Check whitelist/reference preflight as well as normalization, without tasks.
    script = script[:script.index("println groovy.json")] + r'''
 def support = new GroovyClassLoader().parseClass(new File(args[4]))
 support.validateReferenceContract(contract.references, contract.modalities, contract.samples, args[3].toBoolean())
 println groovy.json.JsonOutput.toJson(contract)
 '''
    return subprocess.run(["java", "-cp", str(jar), "groovy.ui.GroovyMain", "-e", script, "--", str(PARSER), str(write_sheet(root, data)), str(root / "results"), str(preprocess).lower(), str(REPO / "lib/WorkflowSupport.groovy")], text=True, capture_output=True, timeout=90)


@pytest.mark.parametrize("kind", FIXTURES)
def test_parser_reference_free_and_full_mode_requirements(tmp_path, kind):
    data = reference_free(fixture(kind))
    result = parse_contract(tmp_path / "preprocess", data)
    assert result.returncode == 0, result.stdout + result.stderr
    contract = json.loads(result.stdout)
    assert contract["samples"]
    assert all("canonical_chromosomes" not in row for row in contract["samples"])
    assert not contract["references"]["dna_bwa_reference"]
    result = parse_contract(tmp_path / "full", data, preprocess=False)
    assert result.returncode != 0
    assert "references.root" in result.stderr


@pytest.mark.parametrize("mutation,error", [
    ("whitelist", "Required ligation barcode whitelist not found"),
    ("mapping", "collision"),
    ("unknown_index", "unknown oligo index"),
    ("missing_fastq", "FASTQ"),
])
def test_preprocess_parser_keeps_input_validation(tmp_path, mutation, error):
    data = reference_free(fixture("rna"))
    sample = next(iter(data["samples"].values()))
    if mutation == "whitelist":
        data["references"]["ligation_barcode_whitelist"] = str(tmp_path / "missing.txt")
    elif mutation == "mapping":
        sample["groups"] = {"A": {"sb_oligo_indices": [5]}, "B": {"sb_oligo_indices": [5]}}
    elif mutation == "unknown_index":
        sample["groups"] = {"A": {"sb_oligo_indices": [9999]}}
    else:
        sample["rna"]["reads"]["r1"] = [str(tmp_path / "missing.fastq")] * len(sample["rna"]["reads"]["r1"])
    result = parse_contract(tmp_path / "parser", data)
    assert result.returncode != 0
    assert error.lower() in result.stderr.lower(), result.stderr


@pytest.fixture(scope="module")
def preprocess_runtime(tmp_path_factory):
    required = ["python3", "codon", "cutadapt", "trim_galore", "pigz", "fastqc", "multiqc"]
    if not all((RUNTIME / "bin" / name).is_file() for name in required):
        pytest.skip("Preprocessing runtime unavailable; set TRESFLOW_TEST_RUNTIME")
    root = tmp_path_factory.mktemp("preprocess_runtime")
    (root / "bin").mkdir()
    for name in required:
        (root / "bin" / name).symlink_to(RUNTIME / "bin" / name)
    (root / "lib").symlink_to(RUNTIME / "lib", target_is_directory=True)
    assert all(not (root / "bin" / name).exists() for name in ALIGNMENT_TOOLS)
    isolated_path = f"{root}/bin:/usr/local/bin:/usr/bin:/bin"
    assert all(shutil.which(name, path=isolated_path) is None for name in ALIGNMENT_TOOLS)
    return root


def workflow(root, data, preprocess=True, cleanup=False, real=False, extra=()):
    nextflow = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
    if not nextflow:
        pytest.skip("Nextflow unavailable")
    sheet = write_sheet(root, data)
    env = os.environ.copy()
    env.update(NXF_OFFLINE="true", NXF_HOME=env.get("NXF_HOME", str(root / "nxf_home")), PATH=f"{data['runtime']['env_prefix']}/bin:/usr/local/bin:/usr/bin:/bin")
    cmd = [nextflow, "-log", str(root / "nextflow.log"), "run", str(REPO), "-ansi-log", "false", "-work-dir", str(root / "work")]
    if not real:
        cmd += ["-profile", "test"]
    cmd += ["--samplesheet", str(sheet), "--outdir", str(root / "results"), "--max_cpus", "4", "--cleanup_work", str(cleanup).lower(), "--preprocess_only", str(preprocess).lower(), "--publish_split_fastqs", "false" if preprocess else "true"]
    result = subprocess.run(cmd + list(extra), cwd=root, env=env, text=True, capture_output=True, timeout=600)
    return result, root / "results"


def trace(out):
    return list(csv.DictReader((out / "pipeline_info/execution_trace.tsv").open(), delimiter="\t"))


def exports(out):
    return {str(p.relative_to(out)): gzip.decompress(p.read_bytes()) for p in out.glob("*_split_fastqs/*.fastq.gz")}


def assert_preprocess(result, out):
    assert result.returncode == 0, result.stdout + result.stderr
    assert "--publish_split_fastqs false is overridden" in result.stdout
    tasks = trace(out)
    assert all(row["status"] == "COMPLETED" for row in tasks)
    assert not {row["name"].split(":")[-1].split(" (")[0] for row in tasks} & DISABLED
    exported = exports(out)
    assert exported
    expected = set()
    for modality in ("rna", "dna"):
        mapping = out / "pipeline_info/derived_contract" / ("rna_sb_group_map.tsv" if modality == "rna" else "dna_mo_map.tsv")
        if not mapping.exists():
            continue
        for row in csv.DictReader(mapping.open(), delimiter="\t"):
            prefix = f"{row['sample']}_{row['sb_group']}" + (f"_{row['mark']}" if modality == "dna" else "")
            expected.update(f"{modality}_split_fastqs/{prefix}_R{mate}.fastq.gz" for mate in (1, 2))
    assert set(exported) == expected
    assert not list(out.glob("*_align"))
    assert not list((out / "pipeline_info/derived_contract").glob("*canonical*"))
    assert "Preprocessing-only" in (out / "TrES_Stats/tres_report.html").read_text()
    qc = list(csv.DictReader((out / "TrES_Stats/qc_metrics.tsv").open(), delimiter="\t"))
    statuses = [r for r in qc if r["metric"].endswith("_status")]
    assert statuses and all(r["value"] == "not run" for r in statuses)
    assert (out / "TrES_Stats/qc/multiqc/multiqc_report.html").is_file()
    assert list((out / "TrES_Stats/qc/fastqc").glob("*.html"))
    for tool in ALIGNMENT_TOOLS:
        assert tool not in {line.split("\t")[0] for line in (out / "pipeline_info/runtime_contract.tsv").read_text().splitlines()}


def assert_pairing(exported):
    for name, raw in exported.items():
        if "_R1." not in name:
            continue
        r1, r2 = raw.decode().splitlines(), exported[name.replace("_R1.", "_R2.")].decode().splitlines()
        assert len(r1) == len(r2) and len(r1) % 4 == 0
        for offset in range(0, len(r1), 4):
            assert r1[offset] == r2[offset]
            assert "NoMatch" not in r1[offset]
            for tag in ("CB:Z:", "XI:Z:", "SB:Z:", "RG:Z:"):
                assert tag in r1[offset]
            assert ("UM:Z:" if name.startswith("rna") else "MO:Z:") in r1[offset]
            tags = dict(token.split(":Z:", 1) for token in r1[offset].split()[1:] if ":Z:" in token)
            assert tags["CB"] == tags["XI"]
            assert re.search(r"_[0-9]+_[ACGT]{24}$", tags["CB"])
            if name.startswith("rna_"):
                primary = r1[offset].split()[0]
                assert primary.endswith("_" + tags["UM"])
                assert re.search(r"_P[0-9a-f]{64}_", primary)
                assert len(r1[offset + 1]) >= 20 and len(r2[offset + 1]) >= 20
            for lines in (r1, r2):
                assert lines[offset + 2] == "+"
                assert len(lines[offset + 1]) == len(lines[offset + 3])


@pytest.mark.parametrize("kind", FIXTURES)
def test_workflow_bulk_exports_preserve_full_mode_and_dna(tmp_path, preprocess_runtime, kind):
    full = extend_rna_fixture(tmp_path, fixture(kind))
    before = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for sample in full["samples"].values() for mod in ("rna", "dna") if mod in sample for paths in sample[mod]["reads"].values() for p in paths}
    # Full-mode baseline includes references and alignment tools.
    result, full_out = workflow(tmp_path / "baseline", full, preprocess=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert {r["name"].split(":")[-1].split(" (")[0] for r in trace(full_out)} & DISABLED
    free = reference_free(full)
    free["runtime"]["env_prefix"] = str(preprocess_runtime)
    result, out = workflow(tmp_path / "preprocess", free)
    assert_preprocess(result, out)
    exported, original = exports(out), exports(full_out)
    assert {k: v for k, v in exported.items() if k.startswith("dna_")} == {k: v for k, v in original.items() if k.startswith("dna_")}
    assert_rna_conversion(exported, original)
    for name, raw in original.items():
        if name.startswith("rna_"):
            for header in raw.decode().splitlines()[::4]:
                assert not re.search(r"_P[0-9a-f]{64}_[ACGTN]+$", header.split()[0])
    filter_tasks = [row for row in trace(out) if "DUAL_TAG_ARTIFACT_FILTER" in row["name"]]
    assert len(filter_tasks) == (1 if kind in {"dual", "combined"} else 0)
    assert_pairing(exports(out))
    result, clean_out = workflow(tmp_path / "cleanup", free, cleanup=True)
    assert_preprocess(result, clean_out)
    assert exports(clean_out) == exports(out)
    assert any("CLEANUP_" in row["name"] for row in trace(clean_out))
    assert not any("CLEANUP_" in row["name"] for row in trace(out))
    assert list((tmp_path / "preprocess/work").rglob("*.fastq"))
    assert not list((tmp_path / "cleanup/work").rglob("*.fastq"))
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest for p, digest in before.items())


def test_zero_read_routing_and_malformed_fastq(tmp_path, preprocess_runtime):
    data = reference_free(fixture("rna"))
    data["runtime"]["env_prefix"] = str(preprocess_runtime)
    sample = next(iter(data["samples"].values()))
    sample["groups"] = {"Empty": {"sb_oligo_indices": [13]}}
    result, out = workflow(tmp_path / "zero", data, cleanup=True)
    assert_preprocess(result, out)
    assert all(raw == b"" for raw in exports(out).values())
    assert "Group-specific routing branch" in (out / "TrES_Stats/read_retention.tsv").read_text()
    bad = tmp_path / "malformed.fastq"
    bad.write_text("@bad\nACGT\n+\n!\n")
    sample["rna"]["reads"]["r1"][0] = str(bad)
    result, _ = workflow(tmp_path / "malformed", data)
    assert result.returncode != 0
    assert "sequence and quality lengths differ" in result.stdout + result.stderr


@pytest.mark.parametrize("value", ["yes", "1", "garbage"])
def test_preprocess_flag_is_strict_boolean(tmp_path, value):
    result, _ = workflow(tmp_path, fixture("rna"), extra=("--preprocess_only", value))
    assert result.returncode != 0
    assert "Parameter --preprocess_only must be true or false" in result.stdout + result.stderr


def test_real_codon_seq_preprocessing_smoke(tmp_path, preprocess_runtime):
    """No test profile: real tagging, trimming, splitting and QC for both modes."""
    data = reference_free(fixture("combined"))
    data["runtime"]["env_prefix"] = str(preprocess_runtime)
    # Bundled genomic reads are 16 nt and intentionally disappear under real
    # --length 20 trimming. Extend only genomic mates in this smoke fixture;
    # retain the actual SB/UMI prefixes and original ligation index reads.
    for sample in data["samples"].values():
        for modality in ("rna", "dna"):
            if modality not in sample:
                continue
            for role in ("r1", "r2"):
                paths = sample[modality]["reads"][role]
                extended = []
                for index, path in enumerate(paths):
                    lines = Path(path).read_text().splitlines()
                    for offset in range(0, len(lines), 4):
                        lines[offset + 1] += "GCTACGTCAGTCCGATGACCTGAGTCGACTGA" * 2
                        lines[offset + 3] += "F" * 64
                    dest = tmp_path / f"{modality}_{role}_{index}.fastq"
                    dest.write_text("\n".join(lines) + "\n")
                    extended.append(str(dest))
                sample[modality]["reads"][role] = extended
    result, out = workflow(tmp_path / "real", data, real=True)
    assert_preprocess(result, out)
    assert_pairing(exports(out))
    # Codon ran in each real tagging/splitting wrapper; artifacts prove that
    # this is distinct from tests run under -profile test.
    commands = "\n".join(p.read_text() for p in (tmp_path / "real/work").rglob(".command.sh"))
    assert "--mock" not in commands
    assert any(exports(out).values())


def test_empty_source_fastqs_still_fail_input_validation(tmp_path, preprocess_runtime):
    data = reference_free(fixture("combined"))
    data["runtime"]["env_prefix"] = str(preprocess_runtime)
    for sample in data["samples"].values():
        for modality in ("rna", "dna"):
            if modality not in sample:
                continue
            for role in sample[modality]["reads"]:
                path = tmp_path / f"empty_{modality}_{role}.fastq"
                path.touch()
                sample[modality]["reads"][role] = [str(path)]
    result, _ = workflow(tmp_path / "empty", data, cleanup=True)
    assert result.returncode != 0
    assert "Synchronized FASTQ streams contain no records" in result.stdout + result.stderr
    assert all(path.is_file() and path.stat().st_size == 0 for path in tmp_path.glob("empty_*.fastq"))


def test_real_trimming_zero_read_branches(tmp_path, preprocess_runtime):
    # Valid raw inputs may lose every pair at the normal --length 20 gate.
    data = reference_free(fixture("combined"))
    data["runtime"]["env_prefix"] = str(preprocess_runtime)
    result, out = workflow(tmp_path / "trimmed_zero", data, real=True, cleanup=True)
    assert_preprocess(result, out)
    assert all(value == b"" for value in exports(out).values())
    assert not list((tmp_path / "trimmed_zero/work").rglob("*.fastq"))
    rows = list(csv.DictReader((out / "TrES_Stats/read_retention.tsv").open(), delimiter="\t"))
    routed = [row for row in rows if row["stage"].endswith("routing branch")]
    assert routed and all(int(row["output_pairs"]) == 0 for row in routed)


def test_rna_export_failure_retains_split_inputs(tmp_path, preprocess_runtime):
    broken = tmp_path / "runtime"
    shutil.copytree(preprocess_runtime, broken, symlinks=True)
    (broken / "bin/pigz").unlink()
    # Tag-record compression also uses pigz. Fail only the export's stdin
    # compression (-c), after the unchanged upstream consumers have succeeded.
    pigz = broken / "bin/pigz"
    real_pigz = shlex.quote(str((preprocess_runtime / "bin/pigz").resolve()))
    pigz.write_text('#!/bin/sh\nif [ "$1" = "-c" ]; then exit 1; fi\nexec ' + real_pigz + ' "$@"\n')
    pigz.chmod(0o755)
    data = reference_free(extend_rna_fixture(tmp_path, fixture("rna")))
    data["runtime"]["env_prefix"] = str(broken)
    result, out = workflow(tmp_path / "failure", data, cleanup=True)
    assert result.returncode != 0
    assert "EXPORT_RNA_FASTQS" in result.stdout + result.stderr
    assert not exports(out)
    split_dirs = [p.parent for p in (tmp_path / "failure/work").rglob(".command.sh") if "run_split_reads_rna.py" in p.read_text()]
    assert split_dirs and all(list(p.glob("*_R1.fastq")) and list(p.glob("*_R2.fastq")) for p in split_dirs)


@pytest.mark.parametrize('configured,max_cpus,expected', [(None, 20, 16), (32, 4, 4), (2, 4, 2)])
def test_rna_export_allocation_is_configurable_and_capped(tmp_path, preprocess_runtime, configured, max_cpus, expected):
    data = reference_free(extend_rna_fixture(tmp_path, fixture('rna')))
    data['runtime']['env_prefix'] = str(preprocess_runtime)
    extra = ['--max_cpus', str(max_cpus)]
    if configured is not None:
        extra += ['--rna_export_cpus', str(configured)]
    result, out = workflow(tmp_path / 'allocation', data, cleanup=False, extra=extra)
    assert_preprocess(result, out)
    export_tasks = [p for p in (tmp_path / 'allocation/work').rglob('.command.sh') if 'export_rna_fastqs.py' in p.read_text()]
    assert len(export_tasks) == 1
    assert f'--threads {expected}' in export_tasks[0].read_text()
    assert f'RNA export CPU budget={expected}:' in export_tasks[0].with_name('.command.err').read_text()
    assert_pairing(exports(out))
