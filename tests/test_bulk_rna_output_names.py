"""Bulk RNA artifact names and read-only compatibility for completed runs."""
import fnmatch
import json
import re
import sys
from pathlib import Path
from unittest import mock

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
from bulk_report import RNA_MATRIX, render, rna_artifact, write_outputs
from bulk_rna import summarize


def test_rna_published_files_are_modality_explicit():
    module = (REPO / 'modules/local/bulk_rna/main.nf').read_text()
    process = module.split('process BULK_RNA_MOLECULES {', 1)[1]
    outputs = [name.replace('${splitName}', 'sample_group')
               for name in re.findall(r'path\("([^"\n]+)"\)', process.split('    script:', 1)[0])]
    pattern = re.search(r'pattern: "([^\"]+)"', process)[1]
    prefix, choices = pattern.split('{', 1)
    patterns = [prefix + suffix for suffix in choices.rstrip('}').split(',')]
    published = {name for name in outputs if any(fnmatch.fnmatch(name, p) for p in patterns)}
    assert published == {'sample_group.rna.dedup.bam', 'sample_group.rna.dedup.bam.bai',
                         'sample_group.rna.rpm.bw', 'sample_group.rna.qc.json'}
    assert 'sample_group.rna.counts.tsv' in outputs
    report_module = (REPO / 'modules/local/bulk_report/main.nf').read_text()
    assert f"path('{RNA_MATRIX}')" in report_module


def test_molecule_summary_writes_explicit_counts_and_preserves_qc_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    index = tmp_path / 'index'
    index.mkdir()
    (index / 'chrName.txt').write_text('chr1\n')
    (index / 'chrStart.txt').write_text('0\n')
    (index / 'exonGeTrInfo.tab').write_text('1\n0\t999\t1\t0\t0\n')
    solo = tmp_path / 'solo'
    raw = solo / 'raw'
    raw.mkdir(parents=True)
    (raw / 'features.tsv').write_text('g1\tgene1\n')
    (raw / 'barcodes.tsv').write_text('pool1\n')
    (raw / 'matrix.mtx').write_text('%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n')
    (solo / 'CellReads.stats').write_text('featureU\tcountedU\n1\t1\n')
    prep = tmp_path / 'preparation.json'
    prep.write_text(json.dumps(dict(routed_pairs=1, retained_pairs=1)))
    (tmp_path / 's_g.Log.final.out').write_text('Uniquely mapped reads number |\t1\nNumber of reads mapped to multiple loci |\t0\n')
    source = tmp_path / 'input.bam'
    with pysam.AlignmentFile(source, 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as bam:
        for flag, pos in [(1 | 64, 10), (1 | 128, 40)]:
            r = pysam.AlignedSegment(bam.header)
            r.query_name = 'pair1'
            r.flag = flag
            r.reference_id = r.next_reference_id = 0
            r.reference_start = pos
            r.next_reference_start = 40 if flag & 64 else 10
            r.query_sequence = 'A' * 10
            r.query_qualities = [30] * 10
            r.cigarstring = '10M'
            r.set_tags([('NH', 1), ('HI', 1), ('GX', 'g1'), ('CB', 'pool1'),
                        ('XP', 'pool1'), ('UB', 'ACGTACGTAC'), ('ZI', 1)])
            bam.write(r)
    summarize('s_g', index, prep, source, solo)
    assert (tmp_path / 's_g.rna.counts.tsv').read_text() == 'gene_id\ts_g\ng1\t1\n'
    assert not (tmp_path / 's_g.counts.tsv').exists()
    assert json.loads((tmp_path / 's_g.rna.qc.json').read_text())['molecules'] == 1


@pytest.mark.parametrize('suffix', ['dedup.bam', 'dedup.bam.bai', 'rpm.bw', 'counts.tsv'])
def test_artifact_lookup_prefers_explicit_names_and_reads_legacy(tmp_path, suffix):
    legacy = tmp_path / f's_g.{suffix}'
    current = tmp_path / f's_g.rna.{suffix}'
    assert rna_artifact(tmp_path, 's_g', suffix) == current
    legacy.write_text('legacy')
    assert rna_artifact(tmp_path, 's_g', suffix) == legacy
    current.write_text('current')
    assert rna_artifact(tmp_path, 's_g', suffix) == current
    assert legacy.read_text() == 'legacy'


@pytest.mark.parametrize('mode', ['new', 'legacy', 'both', 'combined_legacy', 'combined_new', 'combined_both'])
def test_regeneration_preserves_inputs_and_emits_explicit_matrix(tmp_path, mode):
    source = tmp_path / 'inputs'
    source.mkdir()
    expected = 'gene_id\ts_g\ng1\t3\n'
    if mode in {'legacy', 'both'}:
        (source / 's_g.counts.tsv').write_text(expected if mode == 'legacy' else 'gene_id\ts_g\ng1\t99\n')
    if mode in {'new', 'both'}:
        (source / 's_g.rna.counts.tsv').write_text(expected)
    if mode == 'combined_both':
        (source / 'gene_by_group.tsv').write_text('gene_id\ts_g\ng1\t99\n')
    if mode.startswith('combined'):
        (source / ('gene_by_group.tsv' if mode.endswith('legacy') else RNA_MATRIX)).write_text(expected)
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    output = tmp_path / 'regenerated'
    spec = {'targets': [{'id': 's_g', 'modality': 'rna'}]}
    with mock.patch('bulk_report.render', return_value='<html>report</html>'):
        write_outputs(spec, source, output)
    assert (output / RNA_MATRIX).read_text() == expected
    assert not (output / 'gene_by_group.tsv').exists()
    assert {p.name: p.read_bytes() for p in source.iterdir()} == before


def test_report_source_uses_actual_legacy_bam_when_available(tmp_path):
    q = dict(routed_pairs=10, retained_pairs=8, genomeU=6, genomeM=1, featureU=4,
             molecules=2, genes_detected=1, barcodes={tag: {'ACGT': 10} for tag in ('L1', 'L2', 'L3', 'SB')})
    (tmp_path / 's_g.rna.qc.json').write_text(json.dumps(q))
    spec = dict(restart=True, targets=[dict(sample='s', group='g', modality='rna', id='s_g', library='test')])
    assert 's_g.rna.dedup.bam: unavailable' in render(spec, tmp_path)
    (tmp_path / 's_g.dedup.bam').write_text('legacy artifact; not scanned')
    assert 's_g.dedup.bam: unavailable' in render(spec, tmp_path)
    assert not (tmp_path / 's_g.rna.dedup.bam').exists()


def test_dna_only_regeneration_does_not_publish_an_rna_matrix(tmp_path):
    source = tmp_path / 'inputs'
    source.mkdir()
    (source / 'gene_by_group.tsv').write_text('gene_id\ts_g\ng1\t3\n')
    output = tmp_path / 'report'
    with mock.patch('bulk_report.render', return_value='<html>DNA report</html>'):
        write_outputs({'targets': [{'id': 's_g_mark', 'modality': 'dna'}]}, source, output)
    assert not (output / RNA_MATRIX).exists()
    assert (source / 'gene_by_group.tsv').exists()
