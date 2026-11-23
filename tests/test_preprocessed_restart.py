"""Exercise manifest and FASTQ-only restarts with the real toolchain."""
import copy
import csv
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from test_preprocess_only import REPO, RUNTIME, fixture, reference_free, trace, workflow
from test_multi_fastq_inputs import find_nextflow_jar

sys.path.insert(0, str(REPO / "bin"))
from preprocessed_fastqs import validate_split
from tresflow_fastq_utils import fastq_iter


LIGATION = "ACGTACGTTGCATGCAGATCGATC"
REVERSE = str.maketrans("ACGT", "TGCA")


def write_fastq(path, records):
    text = "".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in records)
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def run_tool(name, *args, cwd=None):
    result = subprocess.run([str(RUNTIME / "bin" / name), *map(str, args)],
                            capture_output=True, text=True, cwd=cwd, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def make_fixture(root, intron=False):
    root.mkdir()
    genome = "".join(random.Random(91).choices("ACGT", k=100000))
    if intron:
        genome = genome[:1000] + "GT" + genome[1002:2998] + "AG" + genome[3000:]
    rng = random.Random(43)
    fasta = root / "chr1.fa"
    contigs = {"chr1": genome, **{f"chr{c}": "".join(rng.choices("ACGT", k=1000)) for c in "XYM"}}
    fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in contigs.items()))
    (root / "star").mkdir()
    (root / "bwa").mkdir()
    gtf = root / "genes.gtf"
    exons = [(101, 1000), (3001, 4000)] if intron else [(101, 2000)]
    gtf.write_text("".join(f'chr1\ttest\texon\t{a}\t{b}\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "g1";\n'
                           for a, b in exons))
    run_tool("STAR", "--runMode", "genomeGenerate", "--genomeDir", root / "star",
             "--genomeFastaFiles", fasta, "--sjdbGTFfile", gtf, "--sjdbOverhang", 99,
             "--genomeSAindexNbases", 5, "--genomeChrBinNbits", 10, "--runThreadN", 1, cwd=root)
    run_tool("bwa-mem2", "index", "-p", root / "bwa/chr1", fasta, cwd=root)
    (root / "blacklist.bed").write_text("chr1\t90000\t90001\n")
    (root / "chrom.sizes").write_text("".join(f"{name}\t{len(sequence)}\n" for name, sequence in contigs.items()))
    data = fixture("combined")
    sample = data["samples"]["shared_reads"]
    rna_index = Path(sample["rna"]["reads"]["i1"][0]).read_text().splitlines()[1]
    dna_index = Path(sample["dna"]["reads"]["i1"][0]).read_text().splitlines()[1]
    sample["groups"] = {"shared_pool": {
        "rna_sb_barcodes": ["CAGT"], "dna_sb_barcodes": ["GAT"],
        "mark_barcodes": {"Mark_A": "AGGCTATA"},
    }}
    for modality, index, pairs in (("rna", rna_index, 33), ("dna", dna_index, 48)):
        reads = {}
        for role in ("i1", "r1", "r2"):
            records = []
            for i in range(pairs):
                # Repeated insert coordinates and two physical lanes exercise
                # duplicate marking and optical duplicate isolation.
                name = f"AV240401:AVT0507:TESTFC:{1 + i % 2}:11104:{5031 + i // 2}:3419:ACGT"
                shift = i // 8 * 15
                forward = genome[200 + shift:300 + shift]
                reverse = genome[400 + shift:500 + shift].translate(REVERSE)[::-1]
                if modality == "rna" and i % 4 == 1:
                    # The called-cell BAM needs alignments on both strands to
                    # produce both stranded coverage tracks from this fixture.
                    forward, reverse = reverse, forward
                if role == "i1":
                    seq = index
                elif role == "r1":
                    seq = forward
                else:
                    seq = reverse
                    if modality == "rna":
                        # The final read has enough bases for paired trimming
                        # but fewer than 20 genomic bases after the bulk clip.
                        seq = "ACTG" + "".join(rng.choices("ACGT", k=10)) + seq[:17 if i == 32 else 100]
                records.append((name, seq))
            path = root / f"{modality}_{role}.fastq"
            write_fastq(path, records)
            reads[role] = [str(path)]
        sample[modality]["reads"] = reads
    data["references"].update(
        root=str(root), rna_ref_dir=str(root / "star"), dna_ref_dir=str(root / "bwa"),
        dna_blacklist_bed=str(root / "blacklist.bed"), dna_chrom_sizes=str(root / "chrom.sizes"),
        dna_effective_genome_size=sum(map(len, contigs.values())),
    )
    return data


def records_and_groups(path):
    return (sorted(run_tool("samtools", "view", path).splitlines()),
            sorted(line for line in run_tool("samtools", "view", "-H", path).splitlines()
                   if line.startswith("@RG\t")))


def report_rows(path):
    with path.open() as handle:
        return sorted(tuple(value for field, value in row.items() if field not in {"run", "count_source", "source"})
                      for row in csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize(("options", "exit_code", "message"), [
    (["--help"], 0, "--preprocessed_input true"),
    (["--preprocessed_input", "missing", "--samplesheet", "missing.yaml"], 1, "Restart manifest not found"),
    (["--preprocessed_input", "missing", "--preprocess_only", "true"], 1,
     "Use either --preprocessed_input or --preprocess_only true."),
    (["--preprocessed_input", "true", "--preprocess_only", "true"], 1,
     "Use either --preprocessed_input or --preprocess_only true."),
    (["--restart_from_preprocessed", "missing"], 1, "Use --preprocessed_input"),
])
def test_preprocessed_input_cli(tmp_path, options, exit_code, message):
    nextflow = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
    if not nextflow:
        pytest.skip("Nextflow unavailable")
    environment = dict(os.environ, NXF_OFFLINE="true", NXF_HOME=os.environ.get("NXF_HOME", str(tmp_path / "nxf_home")))
    result = subprocess.run(
        [nextflow, "run", str(REPO), "-ansi-log", "false", *options],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=90,
    )
    output = result.stdout + result.stderr
    assert result.returncode == exit_code and message in output, output
    assert "Submitted process" not in output


def launch_fastq_restart(launch, sheet, opts):
    launch.mkdir()
    nextflow = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
    environment = dict(os.environ, NXF_OFFLINE="true", NXF_HOME=os.environ.get("NXF_HOME", str(launch / "nxf_home")))
    result = subprocess.run(
        [nextflow, "run", str(REPO), "--preprocessed_input", "true", "--samplesheet", str(sheet),
         "--outdir", "results/analysis", "-ansi-log", "false", "--max_cpus", "4",
         "--cleanup_work", "false", *opts],
        cwd=launch, env=environment, capture_output=True, text=True, timeout=600,
    )
    return result, launch / "results/analysis"


def sam_tags(tokens):
    return dict(token.split(":", 2)[::2] for token in tokens)


def assert_rna_restart(r2, out, work, split_id, count):
    expected = {}
    for header, seq, _, qual in fastq_iter(r2):
        name, *comments = header[1:].split()
        expected[name] = (seq, qual, sam_tags(comments))
    assert len(expected) == count

    usam, = [path for path in work.rglob(f"{split_id}_tagged.usam") if not path.is_symlink()]
    sam_records = [line.split("\t") for line in usam.read_text().splitlines() if not line.startswith("@")]
    assert len(sam_records) == 2 * count
    for name, flag, *fields in sam_records:
        seq, qual, tags = expected[name]
        if int(flag) & 128:
            assert fields[7:9] == [seq, qual]
        observed = sam_tags(fields[9:])
        assert all(observed[key] == tags[key] for key in ("CB", "XI", "UM"))

    called = out / f"rna_align/{split_id}.Solo.outGeneFull/raw/barcodes.tsv"
    assert set(called.read_text().splitlines()) == {tags["XI"] for _, _, tags in expected.values()}
    aligned, = [path for path in work.rglob(f"{split_id}.Aligned.sortedByCoord.out.bam") if not path.is_symlink()]
    records = run_tool("samtools", "view", aligned).splitlines()
    assert records
    for line in records:
        name, *fields = line.split("\t")
        tags = sam_tags(fields[10:])
        original = expected[name][2]
        assert tags["CB"] == tags["XI"] == original["XI"]
        assert tags["UM"] == tags["UR"] == tags["UB"] == original["UM"]


def verify_fastq_only_restarts(root, references, runtime, pre, rna, dna, opts, direct):
    for mode in ("rna", "dna", "mixed"):
        active = ("rna", "dna") if mode == "mixed" else (mode,)
        sheet_dir = root / f"sheet-{mode}"
        fastqs = sheet_dir / "recently moved"
        fastqs.mkdir(parents=True)
        sheet = sheet_dir / "preprocessed.yaml"
        entries = []
        retained = {}
        for source in (rna, dna):
            modality = source["modality"]
            if modality not in active:
                continue
            entry = {key: source[key] for key in ("group", "modality")}
            if modality == "dna":
                entry["mark"] = source["mark"]
                if mode == "mixed":
                    entry["mark_barcode"] = "AGGCTATA"
            for role in ("r1", "r2"):
                path = pre / source[role]
                if mode == "mixed" and modality == "rna":
                    path = pre / "rna_split_fastqs" / path.name
                dest = fastqs / f"{modality} renamed {'one' if role == 'r1' else 'two'}.fq.gz"
                retained[dest] = path.read_bytes()
                dest.write_bytes(retained[dest])
                entry[role] = str(dest.relative_to(sheet_dir))
            entries.append(entry)
        data = {
            "library_name": rna["library"],
            "runtime": runtime,
            "references": {key: value for key, value in references.items()
                           if key != "dna_chrom_sizes" and
                           (key in ("root", "species") or any(key.startswith(mod + "_") for mod in active))},
            "samples": {rna["sample"]: {"preprocessed": entries}},
        }
        sheet.write_text(yaml.safe_dump(data))
        launch = root / f"fastq-{mode}"
        result, out = launch_fastq_restart(launch, Path("..") / sheet.relative_to(root), opts)
        assert result.returncode == 0, result.stdout + result.stderr
        assert all(path.read_bytes() == original for path, original in retained.items())
        tasks = {row["name"].split(":")[-1].split(" (")[0] for row in trace(out)}
        assert "VALIDATE_PREPROCESSED_SPLIT" in tasks
        assert not any(name.startswith(("FASTQC", "TAG_", "TRIM_", "SPLIT_RNA_READS",
                                        "SPLIT_DNA_READS", "DUAL_TAG", "CLEANUP_")) for name in tasks)
        assert "Historical preprocessing QC is unavailable" in (out / "TrES_Stats/tres_report.html").read_text()
        with (out / "TrES_Stats/barcode_composition.tsv").open() as handle:
            composition = list(csv.DictReader(handle, delimiter="\t"))
        for modality, split in (("rna", rna), ("dna", dna)):
            if modality in active:
                routed = [row for row in composition if row["modality"] == modality.upper()]
                for tag in ("L1", "L2", "L3"):
                    assert sum(int(row["count"]) for row in routed if row["barcode_type"] == tag) == (
                        split["pairs"] - (modality == "rna" and mode == "mixed"))
                assert all(row["contract_level"] == "routed_sequence" for row in routed)
        assert not (out / "restart_manifest.json").exists()
        assert ("FQ_TO_SAM" in tasks) == ("rna" in active)
        assert ("ALIGN_DNA" in tasks) == ("dna" in active)
        with (out / "TrES_Stats/qc_metrics.tsv").open() as handle:
            metrics = list(csv.DictReader(handle, delimiter="\t"))
        assert {row["parent_sample"] for row in metrics} == {rna["sample"]}
        assert {row["run"] for row in metrics} == {rna["group"]}
        if "rna" in active:
            assert_rna_restart(fastqs / "rna renamed two.fq.gz", out, launch / "work", rna["split_id"],
                               rna["pairs"] - (mode == "mixed"))
            assert any(row["metric"] == "estimated_cells" for row in metrics)
        if "dna" in active:
            assert {"GATK4_MARKDUPLICATES", "DEEPTOOLS_BAMCOVERAGE"} <= tasks
            for stage in ("MarkedDup", "NoDup"):
                path = Path("dna_align") / f"{dna['split_id']}_{stage}.bam"
                assert records_and_groups(out / path) == records_and_groups(direct / path)
            assert list((out / "dna_align").glob("*.bw"))
            assert any(row["branch"] == dna["mark"] and row["metric"] == "read_pairs_examined" for row in metrics)

        if mode == "rna":
            affected = fastqs / "rna renamed two.fq.gz"
            r2_lines = gzip.decompress(retained[affected]).decode().splitlines()
            for failure, header, expected in [
                ("missing-tag", re.sub(r"\s+UM:Z:\S+", "", r2_lines[0]), "missing required FASTQ tags: UM"),
                ("mismatched-mate", r2_lines[0].replace("@", "@other_", 1), "mate read names differ"),
            ]:
                affected.write_bytes(gzip.compress(("\n".join([header, *r2_lines[1:]]) + "\n").encode()))
                failed, _ = launch_fastq_restart(root / failure, sheet, opts)
                log = failed.stdout + failed.stderr
                assert failed.returncode != 0
                assert expected in log
                assert str(affected) in log and rna["split_id"] in log
                if failure == "missing-tag":
                    assert "exact tagged RNA split FASTQs" in log


def test_restart_real_tools(tmp_path):
    required = ("STAR", "bwa-mem2", "samtools", "gatk", "bedGraphToBigWig", "bamCoverage",
                "codon", "trim_galore", "multiqc", "fastqc")
    if not shutil.which("nextflow") and not os.environ.get("NEXTFLOW_BIN"):
        pytest.skip("Nextflow unavailable")
    if not all((RUNTIME / "bin" / name).is_file() for name in required):
        pytest.skip("TrESFlow runtime unavailable")

    data = make_fixture(tmp_path / "inputs")
    opts = ("--tagging_cpus", "1", "--trim_cpus", "1", "--artifact_filter_cpus", "1",
            "--rna_starsolo_cpus", "2", "--dna_align_cpus", "2",
            "--tagging_memory", "2 GB", "--rna_starsolo_memory", "4 GB",
            "--dna_align_memory", "4 GB", "--dna_markdup_memory", "2 GB")
    direct_result, direct = workflow(tmp_path / "direct", data, preprocess=False, real=True, extra=opts)
    assert direct_result.returncode == 0, direct_result.stdout + direct_result.stderr

    pre_result, pre = workflow(tmp_path / "preprocess", reference_free(data), cleanup=True, real=True, extra=opts)
    assert pre_result.returncode == 0, pre_result.stdout + pre_result.stderr
    assert not any(path.is_file() for path in (pre.parent / "work").rglob("*"))
    manifest = json.loads((pre / "restart_manifest.json").read_text())
    assert {(entry["modality"], entry["pairs"]) for entry in manifest["splits"]} == {("rna", 33), ("dna", 48)}
    source_bytes = {}
    for entry in manifest["splits"]:
        for role in ("r1", "r2"):
            relative = Path(entry[role])
            source_bytes[relative] = (pre / relative).read_bytes()
            assert gzip.decompress(source_bytes[relative]) == gzip.decompress(
                (direct / f"{entry['modality']}_split_fastqs" / relative.name).read_bytes())
    dna = next(entry for entry in manifest["splits"] if entry["modality"] == "dna")
    relative = Path(dna["rg"])
    source_bytes[relative] = (pre / relative).read_bytes()
    direct_rg = list((direct.parent / "work").rglob(relative.name))
    assert any(path.read_bytes() == source_bytes[relative] for path in direct_rg)
    rna = next(entry for entry in manifest["splits"] if entry["modality"] == "rna")
    with gzip.open(pre / "rna_split_fastqs" / Path(rna["r2"]).name, "rt") as bulk:
        assert sum(1 for _ in bulk) == (rna["pairs"] - 1) * 4

    restart_sheet = copy.deepcopy(data)
    restart_sheet["references"].pop("ligation_barcode_whitelist")
    for sample in restart_sheet["samples"].values():
        for modality in ("rna", "dna"):
            sample[modality].pop("reads")
    raw_inputs = tmp_path / "inputs"
    saved_reads = raw_inputs / "raw-before-restart"
    saved_reads.mkdir()
    for path in raw_inputs.glob("*.fastq"):
        path.rename(saved_reads / path.name)
    restarted_result, restarted = workflow(
        tmp_path / "restart", restart_sheet, preprocess=False, cleanup=True, real=True,
        extra=(*opts, "--preprocessed_input", str(pre)),
    )
    assert restarted_result.returncode == 0, restarted_result.stdout + restarted_result.stderr
    tasks = {row["name"].split(":")[-1].split(" (")[0] for row in trace(restarted)}
    assert {"VALIDATE_PREPROCESSED_SPLIT", "FQ_TO_SAM", "ALIGN_DNA", "GATK4_MARKDUPLICATES"} <= tasks
    assert not any(name.startswith(("FASTQC", "TAG_", "TRIM_", "SPLIT_RNA_READS",
                                    "SPLIT_DNA_READS", "DUAL_TAG", "CLEANUP_")) for name in tasks)
    for name, before in source_bytes.items():
        assert (pre / name).read_bytes() == before

    matrices = list(direct.glob("rna_align/**/*.mtx"))
    assert matrices
    for path in direct.glob("rna_align/**/*.tsv"):
        assert path.read_bytes() == (restarted / path.relative_to(direct)).read_bytes(), path.name
    for path in matrices:
        assert path.read_bytes() == (restarted / path.relative_to(direct)).read_bytes(), path.name
    for directory in ("rna_align", "dna_align"):
        bams = list((direct / directory).glob("*.bam"))
        bigwigs = list((direct / directory).glob("*.bw"))
        assert bams and bigwigs
        for path in bams:
            assert records_and_groups(path) == records_and_groups(restarted / path.relative_to(direct)), path.name
        for path in bigwigs:
            assert path.read_bytes() == (restarted / path.relative_to(direct)).read_bytes(), path.name
    for path in (direct / "TrES_Stats").glob("*.dna_alignment_retention.tsv"):
        assert path.read_bytes() == (restarted / path.relative_to(direct)).read_bytes()
    for name in ("read_retention.tsv", "qc_metrics.tsv", "barcode_composition.tsv", "library_complexity.tsv"):
        relative = Path("TrES_Stats") / name
        assert report_rows(direct / relative) == report_rows(restarted / relative), name

    verify_fastq_only_restarts(tmp_path, data["references"], data["runtime"], pre, rna, dna, opts, direct)


@pytest.fixture
def dna_split(tmp_path):
    sample = "sample_id"
    group = "group_id"
    name = "I:R:F:1:1000:10:20:UMI"
    cell = f"{sample}_{group}_01_{LIGATION}"
    header = f"{name} CB:Z:{cell}\tXI:Z:{cell}\tSB:Z:GCAT\tMO:Z:AGGCTATA\tRG:Z:I:R:F:L1"
    spec = dict(sample=sample, group=group, library="TestLibrary", modality="dna",
                mo_barcode="AGGCTATA", sb_indices={"GCAT": "01"}, ligation_length=24, pairs=1)
    r1, r2 = tmp_path / "renamed_R1.fastq.gz", tmp_path / "renamed_R2.fastq"
    for path in (r1, r2):
        write_fastq(path, [(header, "ACGTTGCA")])
    return r1, r2, spec, tmp_path / "out.tsv"


def test_regenerate_rg_from_renamed_reads(dna_split, tmp_path):
    r1, r2, spec, output_rg = dna_split
    assert validate_split(r1, r2, spec, output_rg=output_rg) == 1
    expected = "@RG\tID:I:R:F:L1\tSM:sample_id\tLB:TestLibrary\tPU:I:R:F:L1\tPL:ELEMENT\tPM:AVITI_500MIO\n"
    assert output_rg.read_text() == expected
    existing_rg = tmp_path / "source.tsv"
    existing_rg.write_text(expected.replace("TestLibrary", "WrongLibrary"))
    with pytest.raises(ValueError, match="SM/LB/PU"):
        validate_split(r1, r2, spec, existing_rg, output_rg)


def test_fastq_only_label_and_optional_mark_validation(dna_split):
    r1, r2, spec, output_rg = dna_split
    spec = dict(spec, fastq_only=True, sb_indices={}, mo_barcode=None)
    assert validate_split(r1, r2, spec, output_rg=output_rg) == 1
    with pytest.raises(ValueError, match="MO 'AGGCTATA' does not match the DNA mark barcode 'ACGTACGT'"):
        validate_split(r1, r2, dict(spec, mo_barcode="ACGTACGT"), output_rg=output_rg)
    with pytest.raises(ValueError, match="group 'other_group'"):
        validate_split(r1, r2, dict(spec, group="other_group"), output_rg=output_rg)


@pytest.mark.parametrize(("change", "error"), [
    (lambda s: "", "unequal record counts"),
    (lambda s: s.replace("I:R:F:1:1000:", "I:R:F:1:1001:"), "mate read names differ"),
    (lambda s: s.replace("MO:Z:AGGCTATA", "MO:Z:BAD"), "mate FASTQ tags differ"),
    (lambda s: s.replace("IIIIIIII", "IIII"), "sequence and quality lengths differ"),
])
def test_reject_bad_mates(dna_split, change, error):
    r1, r2, spec, output_rg = dna_split
    r2.write_text(change(r2.read_text()))
    with pytest.raises(ValueError, match=error):
        validate_split(r1, r2, spec, output_rg=output_rg)


def test_reject_bulk_rna(tmp_path):
    r1, r2 = tmp_path / "one.fastq", tmp_path / "two.fastq"
    name = "read_P" + "a" * 64 + "_ACGT"
    for path in (r1, r2):
        write_fastq(path, [(name, "ACGT")])
    with pytest.raises(ValueError, match="RNA bulk export"):
        validate_split(r1, r2, dict(modality="rna", ligation_length=24, sb_indices={}))


@pytest.mark.parametrize("missing", ["CB", "XI", "UM"])
def test_clipped_rna_requires_recoverable_cell_identity_and_umi(tmp_path, missing):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    umi = "CCGGTTAACC"
    cell = f"sample_id_group_id_01_{LIGATION}"
    name = f"read_P{'a' * 64}_{umi}"
    tags = {"CB": cell, "XI": cell, "UM": umi}
    for path in (r1, r2):
        comment = "\t".join(f"{key}:Z:{value}" for key, value in tags.items() if key != missing)
        write_fastq(path, [(f"{name} {comment}", "TGCA" * 20)])
    spec = dict(modality="rna", sample="sample_id", group="group_id", fastq_only=True,
                ligation_length=24, umi_start=4, umi_length=10)
    with pytest.raises(ValueError, match=f"missing required FASTQ tags: {missing}.*exact tagged RNA split FASTQs"):
        validate_split(r1, r2, spec)


def test_reject_clipped_rna_without_bulk_name(tmp_path):
    r1, r2 = tmp_path / "one.fastq", tmp_path / "two.fastq"
    cell = f"sample_with_underscores_group_01_{LIGATION}"
    header = f"read CB:Z:{cell}\tXI:Z:{cell}\tRG:Z:{cell}\tSB:Z:GCAT\tUM:Z:CCGGTTAACC"
    spec = dict(modality="rna", sample="sample_with_underscores", group="group",
                ligation_length=24, sb_indices={"GCAT": "01"}, umi_start=4, umi_length=10)
    write_fastq(r1, [(header, "ACGT" * 5)])
    write_fastq(r2, [(header, "GCATGGTTAACCGG" + "ACGT" * 5)])
    assert validate_split(r1, r2, spec) == 1
    write_fastq(r2, [(header, "ACGT" * 5)])
    with pytest.raises(ValueError, match="R2 prefix does not encode"):
        validate_split(r1, r2, spec)


def test_manifest_metadata_and_source_guards(tmp_path):
    jar = find_nextflow_jar()
    if not jar or not shutil.which("java"):
        pytest.skip("Java/local Nextflow jar unavailable")
    for name in ("first.fastq", "second.fastq"):
        (tmp_path / name).touch()
    entry = dict(
        split_id="sample_id_group_id", sample="sample_id", group="group_id",
        modality="rna", library="LIB", format="tresflow-split-v1",
        r1="first.fastq", r2="second.fastq",
    )
    source = tmp_path / "restart_manifest.json"
    source.write_text(json.dumps(dict(version=1, splits=[entry])))
    result = subprocess.run([
        "java", "-cp", str(jar), "groovy.ui.GroovyMain", "-e", """
        def support = new GroovyClassLoader().parseClass(new File(args[0]))
        def base = args[1]
        def rows = [[id: 'sample_id', modality: 'rna', library_name: 'LIB',
                     split_targets: [sample_id_group_id: [group: 'group_id']]]]
        def checkFailure = { String message, Closure operation ->
            try {
                operation()
                assert false: "Expected rejection: " + message
            } catch (IllegalArgumentException e) {
                assert e.message.contains(message): e.message
            }
        }
        def loaded = support.read(base, base + '/out', base + '/work')
        support.bindSamples(loaded, rows)
        assert loaded.splits[0].row.id == 'sample_id'
        checkFailure('outside the restart inputs') { support.read(base, base, base + '/work') }
        checkFailure('outside the restart inputs') { support.read(base, base + '/out', base) }
        loaded.splits[0].entry.group = 'id'
        checkFailure("expected group 'group_id'") { support.bindSamples(loaded, rows) }
        def manifest = new File(base, 'restart_manifest.json')
        manifest.text = manifest.text.replace('"version": 1', '"version": 2')
        checkFailure('requires version 1') { support.read(manifest.path, base + '/out', base + '/work') }
        """, "--", str(REPO / "lib/PreprocessedManifest.groovy"), str(tmp_path)],
        capture_output=True, text=True, timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
