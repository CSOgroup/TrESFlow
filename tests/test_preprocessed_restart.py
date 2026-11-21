"""Validate exact restart inputs and compare a bounded real-toolchain run."""
import copy
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from test_preprocess_only import REPO, RUNTIME, fixture, reference_free, trace, workflow
from test_multi_fastq_inputs import find_nextflow_jar

sys.path.insert(0, str(REPO / "bin"))
from preprocessed_fastqs import validate_split


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


def make_fixture(root):
    root.mkdir()
    genome = "".join(random.Random(91).choices("ACGT", k=100000))
    rng = random.Random(43)
    fasta = root / "chr1.fa"
    contigs = {"chr1": genome, **{f"chr{c}": "".join(rng.choices("ACGT", k=1000)) for c in "XYM"}}
    fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in contigs.items()))
    (root / "star").mkdir()
    (root / "bwa").mkdir()
    gtf = root / "genes.gtf"
    gtf.write_text('chr1\ttest\texon\t101\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "g1";\n')
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
    (["--help"], 0, "--preprocessed_input <directory|manifest.json>"),
    (["--preprocessed_input", "missing", "--samplesheet", "missing.yaml"], 1, "Restart manifest not found"),
    (["--preprocessed_input", "missing", "--preprocess_only", "true"], 1,
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
