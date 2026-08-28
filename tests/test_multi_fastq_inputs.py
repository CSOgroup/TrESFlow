import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
PARSER = REPO / "lib" / "SamplesheetParser.groovy"
TAG_WRAPPER = REPO / "bin" / "run_tag.py"


GROOVY_PARSE = r'''
def parser = new GroovyClassLoader().parseClass(new File(args[0]))
def defaults = [
    rna: [
        sample: [bc_len: 4, bc_start: 0, hd: 0, tag: 'SB', first_pass: 'first_pass', reverse_complement: true],
        umi: [bc_len: 10, bc_start: 4, tag: 'UM'],
        cell: [bc_len: 8, hd: 1, tag: 'CB'],
    ],
    dna: [
        sample: [bc_len: 4, bc_start: 14, hd: 0, tag: 'SB', first_pass: 'first_pass', reverse_complement: true],
        modality: [bc_len: 8, bc_start: 18, hd: 1, tag: 'MO', first_pass: 'not_first_pass', reverse_complement: true],
        cell: [bc_len: 8, hd: 1, tag: 'CB'],
    ],
]
def contract = parser.parseContract(args[1], [outdir: args[2], barcode_defaults: defaults])
println groovy.json.JsonOutput.toJson(contract.samples.collect { row -> [
    id: row.id,
    modality: row.modality,
    i1: row.i1,
    i2: row.i2,
    i2_implicit: row.i2_implicit,
    r1: row.r1,
    r2: row.r2,
] })
'''


def find_nextflow_jar():
    candidates = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.extend(Path(conda_prefix).glob("share/nextflow/dist/*/nextflow-*-one.jar"))
    candidates.extend(Path("/home/annan/micromamba/envs/tres/share/nextflow/dist").glob("*/nextflow-*-one.jar"))
    candidates.extend(Path.home().glob(".nextflow/framework/*/nextflow-*-one.jar"))
    return sorted(candidates)[-1] if candidates else None


def parser_command(sheet, outdir):
    java = shutil.which("java")
    jar = find_nextflow_jar()
    if not java or jar is None:
        pytest.skip("Java and a local Nextflow distribution are required for samplesheet parser tests")
    return [
        java,
        "-cp",
        str(jar),
        "groovy.ui.GroovyMain",
        "-e",
        GROOVY_PARSE,
        "--",
        str(PARSER),
        str(sheet),
        str(outdir),
    ]


def parse_sheet(sheet, outdir, check=True):
    return subprocess.run(
        parser_command(sheet, outdir),
        cwd=REPO,
        text=True,
        capture_output=True,
        check=check,
    )


def fastq_record(name, sequence="ACGT"):
    return f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"


def write_fastq(path, names=("read1",), compressed=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(fastq_record(name) for name in names)
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def base_contract(samples):
    testdata = REPO / "assets" / "testdata"
    return {
        "library_name": "MULTI_TEST",
        "runtime": {"env_prefix": "/tmp/tresflow-test-env", "tmpdir": "/tmp"},
        "references": {
            "species": "human",
            "root": str(testdata / "TrESFlow_References"),
            "ligation_barcode_whitelist": str(
                testdata / "TrESFlow_References" / "ligation_barcode_whitelist.txt"
            ),
            "rna_ref_dir": str(testdata / "TrESFlow_References" / "rna" / "human" / "star"),
            "dna_ref_dir": str(testdata / "TrESFlow_References" / "dna" / "human" / "bwa"),
            "dna_blacklist_bed": str(testdata / "TrESFlow_References" / "dna" / "human" / "blacklist.bed"),
            "dna_chrom_sizes": str(testdata / "TrESFlow_References" / "dna" / "human" / "chrom.sizes"),
            "dna_effective_genome_size": 12,
        },
        "samples": samples,
    }


def rna_sample(reads):
    return {
        "groups": {"group": {"rna_sb_barcodes": ["CAGT"]}},
        "rna": {"reads": reads},
    }


def dna_sample(reads, tagmentation="dual"):
    group = {
        "mark_barcodes": {"mark": "AGGCTATA"},
        "dna_sb_barcodes": ["AAA"],
    }
    if tagmentation == "single":
        group = {"mark_barcodes": {"mark": "TTTACGTA"}, "sb_barcodes": ["CGTA"]}
    return {
        "groups": {"group": group},
        "dna": {"tagmentation": tagmentation, "reads": reads},
    }


def write_sheet(path, samples):
    path.write_text(json.dumps(base_contract(samples), indent=2), encoding="utf-8")


def yaml_quote(value):
    return json.dumps(str(value))


def render_yaml_read_value(paths, syntax):
    quoted = [yaml_quote(path) for path in paths]
    if syntax == "singleton":
        assert len(quoted) == 1
        return quoted[0]
    if syntax == "comma":
        return yaml_quote(",".join(str(path) for path in paths))
    if syntax == "comma_whitespace":
        return yaml_quote("  " + " ,  ".join(str(path) for path in paths) + "  ")
    if syntax == "comma_empty":
        assert len(paths) == 2
        return yaml_quote(f"{paths[0]}, ,{paths[1]}")
    if syntax == "inline":
        return f"[{', '.join(quoted)}]"
    if syntax == "block":
        return "\n" + "\n".join(f"          - {entry}" for entry in quoted)
    raise AssertionError(f"unsupported test YAML syntax: {syntax}")


def write_rna_yaml_sheet(path, reads, syntax_by_role):
    references = base_contract({})["references"]
    read_lines = "\n".join(
        f"        {role}: {render_yaml_read_value(reads[role], syntax_by_role[role])}"
        for role in ("i1", "r1", "r2")
    )
    path.write_text(
        f"""library_name: MULTI_TEST
runtime:
  env_prefix: /tmp/tresflow-test-env
  tmpdir: /tmp
references:
  species: human
  root: {yaml_quote(references['root'])}
  ligation_barcode_whitelist: {yaml_quote(references['ligation_barcode_whitelist'])}
  rna_ref_dir: {yaml_quote(references['rna_ref_dir'])}
  dna_ref_dir: {yaml_quote(references['dna_ref_dir'])}
  dna_blacklist_bed: {yaml_quote(references['dna_blacklist_bed'])}
  dna_chrom_sizes: {yaml_quote(references['dna_chrom_sizes'])}
  dna_effective_genome_size: 12
samples:
  yaml_sample:
    groups:
      group:
        rna_sb_barcodes: [CAGT]
    rna:
      reads:
{read_lines}
""",
        encoding="utf-8",
    )


def normalized_read_sets(row):
    return [
        {role: row[role][index] for role in ("i1", "r1", "r2")}
        for index in range(len(row["i1"]))
    ]


def make_rna_read_paths(root, count=3):
    reads = {role: [] for role in ("i1", "r1", "r2")}
    for role in reads:
        for index in range(count):
            path = root / "reads" / f"{role}.part{index + 1}.fastq.gz"
            write_fastq(path, names=(f"read{index + 1}",), compressed=True)
            reads[role].append(path)
    return reads


def test_yaml_slurper_normalizes_singleton_scalar_to_ordered_lists(tmp_path):
    paths = make_rna_read_paths(tmp_path, count=1)
    relative = {role: [path.relative_to(tmp_path) for path in role_paths] for role, role_paths in paths.items()}
    sheet = tmp_path / "singleton.yaml"
    write_rna_yaml_sheet(sheet, relative, {role: "singleton" for role in relative})

    result = parse_sheet(sheet, tmp_path / "out-singleton")
    row = json.loads(result.stdout.strip().splitlines()[-1])[0]

    assert all(isinstance(row[role], list) for role in ("i1", "r1", "r2"))
    assert normalized_read_sets(row) == [
        {role: str(paths[role][0].resolve()) for role in ("i1", "r1", "r2")}
    ]


def test_yaml_slurper_comma_block_and_inline_lists_have_identical_read_sets(tmp_path):
    paths = make_rna_read_paths(tmp_path, count=3)
    relative = {role: [path.relative_to(tmp_path) for path in role_paths] for role, role_paths in paths.items()}
    representations = {}

    for syntax in ("comma", "comma_whitespace", "block", "inline"):
        sheet = tmp_path / f"{syntax}.yaml"
        write_rna_yaml_sheet(sheet, relative, {role: syntax for role in relative})
        result = parse_sheet(sheet, tmp_path / f"out-{syntax}")
        row = json.loads(result.stdout.strip().splitlines()[-1])[0]
        representations[syntax] = normalized_read_sets(row)

    expected = [
        {role: str(paths[role][index].resolve()) for role in ("i1", "r1", "r2")}
        for index in range(3)
    ]
    assert representations["comma"] == expected
    assert representations["comma_whitespace"] == expected
    assert representations["block"] == expected
    assert representations["inline"] == expected


def test_yaml_slurper_rejects_empty_comma_member(tmp_path):
    paths = make_rna_read_paths(tmp_path, count=2)
    relative = {role: [path.relative_to(tmp_path) for path in role_paths] for role, role_paths in paths.items()}
    sheet = tmp_path / "empty-comma.yaml"
    syntax = {"i1": "comma_empty", "r1": "comma", "r2": "comma"}
    write_rna_yaml_sheet(sheet, relative, syntax)

    result = parse_sheet(sheet, tmp_path / "out-empty", check=False)

    assert result.returncode != 0
    assert "samples.yaml_sample.rna.reads.i1 contains an empty FASTQ path at entry 2" in result.stderr


def test_yaml_slurper_reports_the_missing_list_member(tmp_path):
    paths = make_rna_read_paths(tmp_path, count=3)
    relative = {role: [path.relative_to(tmp_path) for path in role_paths] for role, role_paths in paths.items()}
    relative["i1"][1] = Path("reads/missing.part2.fastq.gz")
    sheet = tmp_path / "missing-list-member.yaml"
    write_rna_yaml_sheet(sheet, relative, {role: "block" for role in relative})

    result = parse_sheet(sheet, tmp_path / "out-missing", check=False)

    assert result.returncode != 0
    assert "samples.yaml_sample.rna.reads.i1 entry 2 FASTQ not found" in result.stderr
    assert str(tmp_path / "reads" / "missing.part2.fastq.gz") in result.stderr


def test_yaml_slurper_reports_rna_read_set_cardinality_mismatch(tmp_path):
    paths = make_rna_read_paths(tmp_path, count=3)
    relative = {role: [path.relative_to(tmp_path) for path in role_paths] for role, role_paths in paths.items()}
    relative["r1"] = relative["r1"][:2]
    sheet = tmp_path / "cardinality-mismatch.yaml"
    write_rna_yaml_sheet(sheet, relative, {role: "inline" for role in relative})

    result = parse_sheet(sheet, tmp_path / "out-cardinality", check=False)

    assert result.returncode != 0
    assert "conflicting technical read-set counts: i1=3, r1=2, r2=3" in result.stderr


def test_parser_accepts_legacy_comma_and_yaml_lists_with_path_edge_cases(tmp_path):
    legacy = {}
    for role in ("i1", "r1", "r2"):
        path = tmp_path / "legacy" / f"legacy_{role}.fastq"
        write_fastq(path)
        legacy[role] = str(path.relative_to(tmp_path))

    multi = {}
    for role in ("i1", "r1", "r2"):
        first = tmp_path / f"{role} first dir" / "same.fastq"
        second_name = (
            "same.fastq"
            if role == "i1"
            else "same,chunk.fastq.gz" if role == "r2" else "same.fastq.gz"
        )
        second = tmp_path / f"{role} second dir" / second_name
        write_fastq(first, names=("read1",))
        write_fastq(second, names=("read2",), compressed=second.suffix == ".gz")
        multi[role] = [str(first), str(second)]

    # Comma scalars trim entries; YAML lists preserve a comma inside r2's filename.
    multi["i1"] = f" {multi['i1'][0]} , {multi['i1'][1]} "
    sheet = tmp_path / "input contract.yaml"
    write_sheet(sheet, {"legacy": rna_sample(legacy), "multi": rna_sample(multi)})
    outdir = tmp_path / "out"
    result = parse_sheet(sheet, outdir)
    rows = json.loads(result.stdout.strip().splitlines()[-1])

    legacy_row = next(row for row in rows if row["id"] == "legacy")
    multi_row = next(row for row in rows if row["id"] == "multi")
    assert all(len(legacy_row[role]) == 1 for role in ("i1", "r1", "r2"))
    assert all(len(multi_row[role]) == 2 for role in ("i1", "r1", "r2"))
    assert Path(multi_row["i1"][0]).name == Path(multi_row["i1"][1]).name == "same.fastq"
    assert multi_row["r2"][1].endswith("same,chunk.fastq.gz")

    provenance = list(
        csv.DictReader(
            (outdir / "pipeline_info" / "derived_contract" / "input_fastq_provenance.tsv").open(),
            delimiter="\t",
        )
    )
    assert len(provenance) == 9
    assert [row["read_set_index"] for row in provenance if row["sample"] == "multi" and row["read_role"] == "i1"] == ["1", "2"]
    assert all(Path(row["canonical_input_path"]).is_absolute() for row in provenance)


def test_dual_dna_implicit_i2_is_positional_and_marked_in_provenance(tmp_path):
    reads = {}
    for role in ("i1", "r1", "r2"):
        paths = []
        for index in range(2):
            path = tmp_path / role / f"chunk{index + 1}.fastq"
            write_fastq(path)
            paths.append(str(path))
        reads[role] = paths
    sheet = tmp_path / "dual.yaml"
    write_sheet(sheet, {"dual": dna_sample(reads)})
    outdir = tmp_path / "out"
    result = parse_sheet(sheet, outdir)
    row = json.loads(result.stdout.strip().splitlines()[-1])[0]
    assert row["i2"] == row["i1"]
    assert row["i2_implicit"] is True
    provenance = (outdir / "pipeline_info" / "derived_contract" / "input_fastq_provenance.tsv").read_text()
    assert provenance.count("implicit_i1_fallback") == 2


def test_dual_dna_explicit_i2_must_match_all_role_counts(tmp_path):
    reads = {}
    for role in ("i1", "i2", "r1", "r2"):
        paths = []
        for index in range(2):
            path = tmp_path / role / f"chunk{index + 1}.fastq"
            write_fastq(path)
            paths.append(str(path))
        reads[role] = paths
    reads["i2"] = reads["i2"][:1]
    sheet = tmp_path / "dual_explicit_i2.yaml"
    write_sheet(sheet, {"dual": dna_sample(reads)})
    result = parse_sheet(sheet, tmp_path / "out", check=False)
    assert result.returncode != 0
    assert "conflicting technical read-set counts: i1=2, i2=1, r1=2, r2=2" in result.stderr


@pytest.mark.parametrize(
    ("mutation", "expected"),
    [
        ("length", "conflicting technical read-set counts: i1=2, r1=1, r2=2"),
        ("empty", "contains an empty FASTQ path at entry 2"),
        ("duplicate", "contains duplicate canonical FASTQ paths"),
        ("suffix", "has an invalid FASTQ suffix"),
        ("directory", "is not a regular file"),
        ("cross_role", "reuses one physical FASTQ in technical read set 1 for roles i1 and r1"),
        ("control", "contains a control character"),
    ],
)
def test_parser_rejects_invalid_read_collections(tmp_path, mutation, expected):
    paths = {}
    for role in ("i1", "r1", "r2"):
        paths[role] = []
        for index in range(2):
            path = tmp_path / role / f"chunk{index + 1}.fastq"
            write_fastq(path)
            paths[role].append(str(path))

    if mutation == "length":
        paths["r1"] = paths["r1"][:1]
    elif mutation == "empty":
        paths["i1"] = f"{paths['i1'][0]}, ,{paths['i1'][1]}"
    elif mutation == "duplicate":
        paths["i1"] = [paths["i1"][0], paths["i1"][0]]
    elif mutation == "suffix":
        invalid = tmp_path / "i1" / "invalid.txt"
        invalid.write_text("not fastq", encoding="utf-8")
        paths["i1"][0] = str(invalid)
    elif mutation == "directory":
        invalid = tmp_path / "directory.fastq"
        invalid.mkdir()
        paths["i1"][0] = str(invalid)
    elif mutation == "cross_role":
        paths["r1"][0] = paths["i1"][0]
    elif mutation == "control":
        paths["i1"][0] = paths["i1"][0] + "\t"

    sheet = tmp_path / f"{mutation}.yaml"
    write_sheet(sheet, {"invalid": rna_sample(paths)})
    result = parse_sheet(sheet, tmp_path / "out", check=False)
    assert result.returncode != 0
    assert expected in result.stderr


def write_manifest(path, entries):
    path.write_text("".join(f"{entry}\n" for entry in entries), encoding="utf-8")


def tag_command(tmp_path, manifests):
    whitelist = tmp_path / "whitelist.txt"
    whitelist.write_text("AAAA\n", encoding="utf-8")
    return [
        sys.executable,
        str(TAG_WRAPPER),
        "--mode",
        "mock",
        "--script",
        str(REPO / "scripts" / "core_runtime" / "Tag.codon"),
        "--i2-manifest",
        str(manifests["i2"]),
        "--r1-manifest",
        str(manifests["r1"]),
        "--r2-manifest",
        str(manifests["r2"]),
        "--whitelist",
        str(whitelist),
        "--sample",
        "sample",
        "--tag",
        "SB",
        "--bc-len",
        "4",
        "--bc-start",
        "0",
        "--hd",
        "0",
        "--first-pass-arg",
        "first_pass",
        "--rev-comp-arg",
        "fw",
        "--output-r1",
        str(tmp_path / "tagged.R1.fastq"),
        "--output-r2",
        str(tmp_path / "tagged.R2.fastq"),
        "--output-counts",
        str(tmp_path / "counts.tsv"),
        "--output-stats",
        str(tmp_path / "stats.tsv"),
    ]


def make_virtual_streams(tmp_path, role_records):
    manifests = {}
    for role, chunks in role_records.items():
        entries = []
        for index, records in enumerate(chunks):
            directory = tmp_path / f"{role} source {index + 1},dir"
            suffix = ".fastq.gz" if index % 2 else ".fastq"
            path = directory / f"identical{suffix}"
            write_fastq(path, names=records, compressed=path.suffix == ".gz")
            entries.append(path)
        manifest = tmp_path / f"{role}.manifest"
        write_manifest(manifest, entries)
        manifests[role] = manifest
    return manifests


def test_mock_tagger_streams_mixed_compression_in_manifest_order(tmp_path):
    manifests = make_virtual_streams(
        tmp_path,
        {
            "i2": [("read1",), ("read2",)],
            "r1": [("read1/1",), ("read2/1",)],
            "r2": [("read1/2",), ("read2/2",)],
        },
    )
    # Replace barcode reads with accepted AAAA sequences while preserving QNAMEs.
    for index, path in enumerate((tmp_path / "i2.manifest").read_text().splitlines()):
        write_fastq(Path(path), names=(f"read{index + 1}",), compressed=path.endswith(".gz"))
        if path.endswith(".gz"):
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(fastq_record(f"read{index + 1}", "AAAA"))
        else:
            Path(path).write_text(fastq_record(f"read{index + 1}", "AAAA"), encoding="utf-8")

    subprocess.run(tag_command(tmp_path, manifests), cwd=REPO, check=True, capture_output=True, text=True)
    headers = (tmp_path / "tagged.R1.fastq").read_text().splitlines()[::4]
    stats = (tmp_path / "stats.tsv").read_text()
    assert [header.split()[0] for header in headers] == ["@read1/1", "@read2/1"]
    assert "reads\t2\t100.0%" in stats
    assert (tmp_path / "counts.tsv").read_text().strip() == "2\tAAAA"


def test_mock_tagger_retains_legacy_single_path_cli(tmp_path):
    inputs = {}
    for role in ("i2", "r1", "r2"):
        path = tmp_path / f"{role}.fastq"
        sequence = "AAAA" if role == "i2" else "ACGT"
        path.write_text(fastq_record("legacy", sequence), encoding="utf-8")
        inputs[role] = path
    manifests = make_virtual_streams(
        tmp_path / "unused",
        {role: [("unused",)] for role in ("i2", "r1", "r2")},
    )
    command = tag_command(tmp_path, manifests)
    for role in ("i2", "r1", "r2"):
        manifest_flag = f"--{role}-manifest"
        index = command.index(manifest_flag)
        command[index:index + 2] = [f"--{role}", str(inputs[role])]
    subprocess.run(command, cwd=REPO, check=True, capture_output=True, text=True)
    assert (tmp_path / "tagged.R1.fastq").read_text().startswith("@legacy")


@pytest.mark.parametrize(
    ("case", "expected"),
    [
        ("truncated", "Incomplete FASTQ record"),
        ("unequal", "unequal EOF"),
        ("qname", "QNAME mismatch"),
    ],
)
def test_mock_tagger_rejects_incomplete_or_unsynchronized_streams(tmp_path, case, expected):
    role_records = {role: [("read1",)] for role in ("i2", "r1", "r2")}
    if case == "unequal":
        role_records["i2"] = [("read1", "read2")]
    elif case == "qname":
        role_records["r2"] = [("different",)]
    manifests = make_virtual_streams(tmp_path, role_records)
    if case == "truncated":
        r2_path = Path(manifests["r2"].read_text().strip())
        r2_path.write_text("@read1\nACGT\n+\n", encoding="utf-8")
    result = subprocess.run(tag_command(tmp_path, manifests), cwd=REPO, capture_output=True, text=True)
    assert result.returncode != 0
    assert expected in result.stderr


def test_mock_tagger_does_not_synchronize_across_read_set_boundaries(tmp_path):
    manifests = make_virtual_streams(
        tmp_path,
        {
            "i2": [("read1",), ("read2", "read3")],
            "r1": [("read1", "read2"), ("read3",)],
            "r2": [("read1", "read2"), ("read3",)],
        },
    )
    result = subprocess.run(tag_command(tmp_path, manifests), cwd=REPO, capture_output=True, text=True)
    assert result.returncode != 0
    assert "unequal EOF in technical read set 1" in result.stderr


def test_production_codon_taggers_accept_ordered_manifests(tmp_path):
    codon = shutil.which("codon")
    if not codon:
        pytest.skip("Codon is required for the production manifest integration test")

    names = ("codon_read1", "codon_read2")
    sources = {}
    for role in ("index", "r1", "r2", "ligation"):
        entries = []
        for index, name in enumerate(names):
            suffix = ".fastq.gz" if index else ".fastq"
            path = tmp_path / f"{role} source {index + 1}" / f"same{suffix}"
            if role == "index":
                sequence = "AAAA" if index == 0 else "CCCC"
                record_name = name
            elif role == "ligation":
                sequence = "ACGTACGTACGT"
                record_name = name
            else:
                sequence = "TTTTACGTACGTACGTACGT"
                record_name = f"{name}/{'1' if role == 'r1' else '2'}"
            path.parent.mkdir(parents=True, exist_ok=True)
            text = fastq_record(record_name, sequence)
            if path.suffix == ".gz":
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.write(text)
            else:
                path.write_text(text, encoding="utf-8")
            entries.append(path)
        manifest = tmp_path / f"{role}.manifest"
        write_manifest(manifest, entries)
        sources[role] = manifest

    env = os.environ.copy()
    env["TMPDIR"] = str(tmp_path)
    tag_whitelist = tmp_path / "tag_whitelist.txt"
    tag_whitelist.write_text("AAAA\nCCCC\n", encoding="utf-8")
    tagged_r1 = tmp_path / "tagged.R1.fastq"
    tagged_r2 = tmp_path / "tagged.R2.fastq"
    read_set_counts = tmp_path / "read_set_counts.tsv"
    tag_result = subprocess.run(
        [
            sys.executable,
            str(TAG_WRAPPER),
            "--mode", "real",
            "--script", str(REPO / "scripts/core_runtime/Tag.codon"),
            "--i2-manifest", str(sources["index"]),
            "--r1-manifest", str(sources["r1"]),
            "--r2-manifest", str(sources["r2"]),
            "--whitelist", str(tag_whitelist),
            "--sample", "codon_multi",
            "--tag", "SB",
            "--bc-len", "4",
            "--bc-start", "0",
            "--hd", "0",
            "--first-pass-arg", "first_pass",
            "--rev-comp-arg", "fw",
            "--output-r1", str(tagged_r1),
            "--output-r2", str(tagged_r2),
            "--output-counts", str(tmp_path / "tag.counts.tsv"),
            "--output-stats", str(tmp_path / "tag.stats.tsv"),
            "--output-read-set-counts", str(read_set_counts),
        ],
        cwd=REPO,
        env=env,
        text=True,
        capture_output=True,
    )
    assert tag_result.returncode == 0, tag_result.stdout + tag_result.stderr

    tagged_r1_manifest = tmp_path / "tagged_r1.manifest"
    tagged_r2_manifest = tmp_path / "tagged_r2.manifest"
    write_manifest(tagged_r1_manifest, [tagged_r1])
    write_manifest(tagged_r2_manifest, [tagged_r2])
    umi_r1 = tmp_path / "umi.R1.fastq"
    umi_r2 = tmp_path / "umi.R2.fastq"
    umi_result = subprocess.run(
        [
            sys.executable,
            str(REPO / "bin/run_tag_umi.py"),
            "--mode", "real",
            "--script", str(REPO / "scripts/core_runtime/Tag_UMI.codon"),
            "--i2-manifest", str(sources["r2"]),
            "--r1-manifest", str(tagged_r1_manifest),
            "--r2-manifest", str(tagged_r2_manifest),
            "--sample", "codon_multi",
            "--tag", "UM",
            "--bc-len", "4",
            "--bc-start", "4",
            "--output-r1", str(umi_r1),
            "--output-r2", str(umi_r2),
            "--output-counts", str(tmp_path / "umi.counts.tsv"),
            "--read-set-counts", str(read_set_counts),
        ],
        cwd=REPO,
        env=env,
        text=True,
        capture_output=True,
    )
    assert umi_result.returncode == 0, umi_result.stdout + umi_result.stderr

    umi_r1_manifest = tmp_path / "umi_r1.manifest"
    umi_r2_manifest = tmp_path / "umi_r2.manifest"
    write_manifest(umi_r1_manifest, [umi_r1])
    write_manifest(umi_r2_manifest, [umi_r2])
    ligation_whitelist = tmp_path / "ligation_whitelist.txt"
    ligation_whitelist.write_text("ACGT\n", encoding="utf-8")
    lig_result = subprocess.run(
        [
            sys.executable,
            str(REPO / "bin/run_tag_lig3.py"),
            "--mode", "real",
            "--script", str(REPO / "scripts/core_runtime/Tag_Lig3.codon"),
            "--i1-manifest", str(sources["ligation"]),
            "--r1-manifest", str(umi_r1_manifest),
            "--r2-manifest", str(umi_r2_manifest),
            "--whitelist", str(ligation_whitelist),
            "--sample", "codon_multi",
            "--tag", "CB",
            "--bc-len", "4",
            "--hd", "0",
            "--start-positions", "0,4,8",
            "--output-r1", str(tmp_path / "cell.R1.fastq"),
            "--output-r2", str(tmp_path / "cell.R2.fastq"),
            "--output-counts", str(tmp_path / "cell.counts.tsv"),
            "--output-tag-records", str(tmp_path / "tag_records.tsv"),
            "--read-set-counts", str(read_set_counts),
            "--output-stats", str(tmp_path / "cell.L1.stats.tsv"),
            "--output-stats", str(tmp_path / "cell.L2.stats.tsv"),
            "--output-stats", str(tmp_path / "cell.L3.stats.tsv"),
        ],
        cwd=REPO,
        env=env,
        text=True,
        capture_output=True,
    )
    assert lig_result.returncode == 0, lig_result.stdout + lig_result.stderr
    assert read_set_counts.read_text().splitlines() == [
        "read_set_index\trecord_count",
        "1\t1",
        "2\t1",
    ]
    assert (tmp_path / "tag.stats.tsv").read_text().startswith("reads\t2\t")
    assert len((tmp_path / "cell.R1.fastq").read_text().splitlines()) == 8
    assert len((tmp_path / "tag_records.tsv").read_text().splitlines()) == 2


def test_production_codon_taggers_retain_legacy_single_fastq_cli(tmp_path):
    codon = shutil.which("codon")
    if not codon:
        pytest.skip("Codon is required for the legacy standalone integration test")

    i2 = tmp_path / "legacy_I2.fastq"
    r1 = tmp_path / "legacy_R1.fastq"
    r2 = tmp_path / "legacy_R2.fastq"
    i2.write_text(fastq_record("legacy", "AAAA"), encoding="utf-8")
    r1.write_text(fastq_record("legacy/1", "TTTTACGTACGTACGTACGT"), encoding="utf-8")
    r2.write_text(fastq_record("legacy/2", "TTTTACGTACGTACGTACGT"), encoding="utf-8")
    tag_whitelist = tmp_path / "tag_whitelist.txt"
    tag_whitelist.write_text("AAAA\n", encoding="utf-8")

    tag_out = tmp_path / "tag_out"
    tag_out.mkdir()
    subprocess.run(
        [
            codon, "run", "-plugin", "seq", "-release",
            "-D", "BC_LEN=4", "-D", "BC_START=0", "-D", "HD=0",
            str(REPO / "scripts/core_runtime/Tag.codon"),
            str(i2), str(r1), str(r2), str(tag_whitelist), "legacy", "SB",
            str(tag_out), "first_pass", "fw",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    tagged_r1 = tag_out / "legacy_R1_SB.fastq"
    tagged_r2 = tag_out / "legacy_R2_SB.fastq"
    assert tagged_r1.exists() and tagged_r2.exists()

    umi_out = tmp_path / "umi_out"
    umi_out.mkdir()
    subprocess.run(
        [
            codon, "run", "-plugin", "seq", "-release",
            "-D", "BC_LEN=4", "-D", "BC_START=4",
            str(REPO / "scripts/core_runtime/Tag_UMI.codon"),
            str(r2), str(tagged_r1), str(tagged_r2), "legacy", "UM", str(umi_out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    umi_r1 = umi_out / "legacy_R1_SB_UM.fastq"
    umi_r2 = umi_out / "legacy_R2_SB_UM.fastq"
    assert umi_r1.exists() and umi_r2.exists()

    ligation_sequence = list("A" * 99)
    for start in (15, 53, 91):
        ligation_sequence[start:start + 4] = "ACGT"
    ligation = tmp_path / "legacy_I1.fastq"
    ligation.write_text(fastq_record("legacy", "".join(ligation_sequence)), encoding="utf-8")
    ligation_whitelist = tmp_path / "ligation_whitelist.txt"
    ligation_whitelist.write_text("ACGT\n", encoding="utf-8")
    lig_out = tmp_path / "lig_out"
    lig_out.mkdir()
    subprocess.run(
        [
            codon, "run", "-plugin", "seq", "-release",
            "-D", "BC_LEN=4", "-D", "HD=0",
            str(REPO / "scripts/core_runtime/Tag_Lig3.codon"),
            str(ligation), str(umi_r1), str(umi_r2), str(ligation_whitelist),
            "legacy", "CB", str(lig_out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert (lig_out / "legacy_R1_SB_UM_CB.fastq").exists()
    assert (lig_out / "legacy_R2_SB_UM_CB.fastq").exists()
