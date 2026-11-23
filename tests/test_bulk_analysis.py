"""Small biological equivalence fixture. All input, reference and work paths are private to the test."""
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pyBigWig
import pysam
import pytest
import yaml

from test_preprocessed_restart import REPO, RUNTIME, REVERSE, make_fixture, write_fastq
from test_preprocess_only import reference_free
from tresflow_fastq_utils import fastq_iter
from export_rna_fastqs import pool_identifier
from prepare_bulk_rna import prepare


def controlled(root):
    sheet = make_fixture(root, intron=True)
    genome = (root / "chr1.fa").read_text().splitlines()[1]
    sample = sheet["samples"]["shared_reads"]
    sample["groups"]["shared_pool"]["rna_sb_barcodes"] += ["ACGT"]
    sample["groups"]["other"] = {"rna_sb_barcodes": ["GCAT"]}
    # The first group has three molecules. U1 and its 1MM neighbour share a pool,
    # U1 occurs in a second physical pool, and U2 aligns entirely within an intron.
    u1, u2 = "ACGTACGTAC", "TGCATGCATG"
    options = [(200, "CAGT", u1), (1500, "CAGT", u2), (220, "ACGT", u1),
               (320, "GCAT", u1), (350, "AAAA", u1), (300, "CAGT", u1),
               (400, "CAGT", u1[:-1] + "G"), (230, "CAGT", u1)]
    index = Path(sample["rna"]["reads"]["i1"][0]).read_text().splitlines()[1]
    for role in ("i1", "r1", "r2"):
        records = []
        for i, (pos, sb, umi) in enumerate(options):
            seq = {"i1": index, "r1": genome[pos:pos+80],
                   "r2": (umi + sb).translate(REVERSE)[::-1] +
                   genome[pos+120:pos+(137 if i == 7 else 200)].translate(REVERSE)[::-1]}[role]
            records.append((f"fragment_{i}", seq))
        paths = []
        for i, chunk in enumerate((records[:5], records[5:])):
            path = root / f"lane_{i}_{role}.fq"
            write_fastq(path, chunk)
            paths.append(str(path))
        sample["rna"]["reads"][role] = paths
    return sheet


def launch(root, args):
    environment = dict(os.environ, NXF_OFFLINE="true", NXF_WORK=str(root / ("work-" + args[-1].split("/")[-1])))
    result = subprocess.run([os.environ.get("NEXTFLOW_BIN", "nextflow"), "run", ".", *args],
                            cwd=root, env=environment, capture_output=True, text=True, timeout=900)
    (root / ("run-" + args[-1].split("/")[-1] + ".log")).write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr


def project(root):
    root.mkdir()
    for item in ("main.nf", "bin", "modules", "subworkflows", "conf", "scripts", "assets", "workflows", "lib"):
        (root / item).symlink_to(REPO / item)
    # Relative includes resolve from this isolated launch project.
    (root / "nextflow.config").write_text((REPO / "nextflow.config").read_text() + """
params.max_cpus = 4
params.tagging_cpus = 1
params.tagging_memory = '2 GB'
params.rna_starsolo_cpus = 2
params.rna_starsolo_memory = '4 GB'
params.dna_align_cpus = 2
params.dna_align_memory = '4 GB'
params.dna_markdup_memory = '2 GB'
""")


def snapshot(root):
    result = {}
    for p in root.rglob("*"):
        if p.suffix == ".bam":
            with pysam.AlignmentFile(p) as bam:
                assert bam.has_index()
                result[p.name] = [r.to_string() for r in bam]
        elif p.suffix == ".bw":
            with pyBigWig.open(str(p)) as bw:
                result[p.name] = {c: bw.intervals(c) for c in bw.chroms()}
        elif p.suffix == ".bed" or p.name == "gene_by_group.tsv":
            result[p.name] = p.read_text()
        elif p.name.endswith(".qc.json"):
            qc = json.loads(p.read_text())
            # The one discarded short RNA pair is already absent from clipped exports.
            for key in ("routed_pairs", "barcodes"):
                qc.pop(key, None)
            result[p.name] = qc
    return result


def test_bulk_routes(tmp_path):
    if not (RUNTIME / "bin/epic2").exists() or not shutil.which("nextflow"):
        pytest.skip("bulk runtime unavailable")
    root = tmp_path / "pipeline"
    project(root)
    raw = controlled(tmp_path / "inputs")
    (root / "raw.yaml").write_text(yaml.safe_dump(raw))
    (root / "preprocess.yaml").write_text(yaml.safe_dump(reference_free(raw)))
    launch(root, ["--bulk", "true", "--samplesheet", "raw.yaml", "--outdir", "results/bulk_raw"])
    start = root / "results/bulk_raw/bulk"
    with (start / "gene_by_group.tsv").open() as h:
        assert list(csv.DictReader(h, delimiter="\t")) == [
            {"gene_id": "g1", "shared_reads_shared_pool": "3", "shared_reads_other": "1"}]
    q = json.loads((start / "rna/shared_reads/shared_pool/shared_reads_shared_pool.rna.qc.json").read_text())
    assert q["intronic"] == 1 and q["sense"] == 5 and q["pools_detected"] == 2
    dna, = start.rglob("*.dna.qc.json")
    d = json.loads(dna.read_text())
    assert d["usable_duplicate_fragments"] == 42 and d["usable_fragments"] == 48
    assert d["frip"] == 1 and d["peak_count"] == 1
    track, = start.rglob("*.cpm.bw")
    with pyBigWig.open(str(track)) as bw:
        assert bw.values("chr1", 200, 201)[0] == pytest.approx(32 / 48 * 1e6, rel=1e-5)
    assert "had no sample-barcode group" in (start / "bulk_report.html").read_text()
    launch(root, ["--preprocess_only", "true", "--samplesheet", "preprocess.yaml", "--outdir", "results/pre"])
    pre = root / "results/pre"
    copies = root / "sheets" / "pairs"
    copies.mkdir(parents=True)
    manifest = json.loads((pre / "restart_manifest.json").read_text())
    for fmt, out in (("exact", "bulk_exact"), ("clipped", "bulk_from_fastqs")):
        entries, hashes = [], {}
        for row in manifest["splits"]:
            item = {k: row[k] for k in ("modality", "group", "mark") if k in row}
            for mate in ("r1", "r2"):
                source = pre / row[mate]
                if fmt == "clipped" and row["modality"] == "rna":
                    source = pre / "rna_split_fastqs" / source.name
                copy = copies / f"{len(entries)} renamed {fmt} {mate}.fq.gz"
                shutil.copyfile(source, copy)
                hashes[copy] = hashlib.sha256(copy.read_bytes()).digest()
                item[mate] = "pairs/" + copy.name
            entries.append(item)
        sheet = dict(raw, samples={"shared_reads": {"preprocessed": entries}}, references=dict(raw["references"]))
        sheet["references"].pop("ligation_barcode_whitelist")
        path = root / "sheets" / "preprocessed.yaml"
        path.write_text(yaml.safe_dump(sheet))
        # The samplesheet may itself be a symlink. Use an actual nested path so
        # resolution of pairs/ is demonstrably independent of the launch folder.
        launch(root, ["--bulk", "true", "--preprocessed_input", "true", "--samplesheet",
                      "sheets/preprocessed.yaml", "--outdir", f"results/{out}"])
        result = root / "results" / out / "bulk"
        assert snapshot(start) == snapshot(result)
        assert "unassigned totals" in (result / "bulk_report.html").read_text()
        assert all(hashlib.sha256(p.read_bytes()).digest() == digest for p, digest in hashes.items())
        assert not any(p.is_file() for p in (root / f"work-{out}").rglob("*"))


def test_shared_logical_identity_uses_physical_sb(tmp_path):
    spec = dict(sample="s", group="g", modality="rna", fastq_only=True, ligation_length=24,
                umi_start=4, umi_length=10, sb_start=0, sb_length=4, sb_hd=0, sb_rc=True)
    umi = "ACGTACGTAC"
    cell = "s_g_01_" + "ACGT" * 6
    records = [[], []]
    for i, sb in enumerate(("ATGC", "CGAT")):
        for mate in (0, 1):
            prefix = (umi + sb).translate(REVERSE)[::-1] if mate else ""
            seq = prefix + "GACT" * 6
            records[mate].append((f"same{i} CB:Z:{cell} XI:Z:{cell} UM:Z:{umi}", seq))
    r1, r2 = (tmp_path / f"mate{n}.fq" for n in (1, 2))
    for p, reads in zip((r1, r2), records):
        write_fastq(p, reads)
    output = tmp_path / "group.sam"
    qc = prepare(spec, r1, r2, output)
    assert qc["retained_pairs"] == 2
    assert {line.split("XP:Z:")[1].split()[0] for line in output.read_text().splitlines()[1:]} == {
        pool_identifier(cell, sb) for sb in ("ATGC", "CGAT")}
    with pytest.raises(ValueError, match="missing physical SB"):
        prepare(dict(spec, sb_hd=1), r1, r2, tmp_path / "ambiguous.sam")


def test_dna_only_single_targets_and_reference(tmp_path):
    root = tmp_path / "pipeline"
    project(root)
    inputs = tmp_path / "inputs"
    sheet = make_fixture(inputs)
    sample = sheet["samples"]["shared_reads"]
    del sample["rna"]
    sample["groups"] = {"pool": {"dna_sb_barcodes": ["CGTA"], "mark_barcodes": {"target": "TTTACGTA"}}}
    dna = sample["dna"]
    dna["tagmentation"] = "single"
    records = list(fastq_iter(Path(dna["reads"]["r1"][0])))
    for role, filename in (("i1", "test_dna_I1.fastq"), ("i2", "test_dna_I2_single_lig.fastq")):
        seq = (REPO / "assets/testdata" / filename).read_text().splitlines()[1]
        path = inputs / (role + ".fq")
        write_fastq(path, [(r[0][1:], seq) for r in records])
        dna["reads"][role] = str(path)
    # Retain all 32 fragments overlapping bin 200 and just one at each of
    # the last two positions. CPM now distinguishes 32/34 from nodup 4/6.
    for role, paths in dna["reads"].items():
        path = Path(paths[0] if isinstance(paths, list) else paths)
        reads = list(fastq_iter(path))
        write_fastq(path, [(r[0][1:], r[1]) for i, r in enumerate(reads) if i < 32 or i in (32, 40)])
    (root / "overlap.bed").write_text("chr1\t200\t250\nchr1\t210\t240\n")
    (root / "raw.yaml").write_text(yaml.safe_dump(sheet))
    launch(root, ["--bulk", "true", "--bulk_reference_bed", "overlap.bed", "--samplesheet", "raw.yaml",
                  "--outdir", "results/single"])
    bulk = root / "results/single/bulk"
    qc = json.loads(next(bulk.rglob("*.dna.qc.json")).read_text())
    assert qc["usable_fragments"] == 34 and qc["usable_duplicate_fragments"] == 28
    assert qc["fragments_in_reference"] == 32
    with pyBigWig.open(str(next(bulk.rglob("*.bw")))) as track:
        assert track.values("chr1", 200, 201)[0] == pytest.approx(32/34 * 1e6, rel=1e-5)
    assert not (bulk / "gene_by_group.tsv").exists()


def test_rna_only_physical_pools_and_splice(tmp_path):
    root = tmp_path / "pipeline"
    project(root)
    inputs = tmp_path / "inputs"
    sheet = controlled(inputs)
    genome = (inputs / "chr1.fa").read_text().splitlines()[1]
    cell, umi = "s_g_01_" + "ACGT" * 6, "ACGTACGTAC"
    prefix = "TGCATGCATG"
    first, second = [], []
    for i, (sb, molecular, r1, r2) in enumerate([
        ("ATGC", umi, genome[960:1000] + genome[3000:3040], genome[3050:3130]),
        ("CGAT", umi, genome[400:480], genome[540:620]),
        # STAR counts this unique molecule from a paired template where only
        # R1 aligns. Both records must survive, with truthful mapping flags.
        ("CGAT", prefix, genome[600:680], "N" * 80),
    ]):
        tags = f"CB:Z:{cell} XI:Z:{cell} UM:Z:{molecular} SB:Z:{sb}"
        name = f"read{i}_{pool_identifier(cell, sb)}_{molecular} {tags}"
        first.append((name, r1))
        second.append((name, r2.translate(REVERSE)[::-1]))
    for role, records in (("r1", first), ("r2", second)):
        write_fastq(root / f"{role}.fq", records)
    sheet["samples"] = {"s": {"preprocessed": [{"modality": "rna", "group": "g", "r1": "r1.fq", "r2": "r2.fq"}]}}
    (root / "preprocessed.yaml").write_text(yaml.safe_dump(sheet))
    launch(root, ["--bulk", "true", "--preprocessed_input", "true", "--samplesheet", "preprocessed.yaml",
                  "--outdir", "results/bulk_from_fastqs"])
    out = root / "results/bulk_from_fastqs/bulk"
    with pysam.AlignmentFile(next(out.rglob("*.bam"))) as bam:
        records = list(bam)
    assert len(records) == 6 and sum(r.is_unmapped for r in records) == 1
    pools = {r.get_tag("CB") for r in records if r.get_tag("UB") == umi}
    assert pools == {pool_identifier(cell, sb) for sb in ("ATGC", "CGAT")}
    assert all(r.get_tag("XI") == cell for r in records)
    assert any("N" in (r.cigarstring or "") for r in records)
    with pyBigWig.open(str(next(out.rglob("*.bw")))) as bw:
        import math
        assert bw.values("chr1", 980, 981)[0] > 0
        gap = bw.values("chr1", 2000, 2001)[0]
        assert math.isnan(gap) or gap == 0
    assert (out / "gene_by_group.tsv").read_text() == "gene_id\ts_g\ng1\t3\n"
