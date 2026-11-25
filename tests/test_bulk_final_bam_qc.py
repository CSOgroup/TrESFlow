"""Audit real final BAM representatives through producer QC and the bulk report."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pysam

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
from bulk_report import render, write_outputs
from bulk_rna import summarize, update_final_bam_qc


def make_record(header, name, flag, position, umi, identity):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.flag = flag
    read.query_sequence = 'A' * 20
    read.query_qualities = [30] * 20
    read.reference_id = read.next_reference_id = 0
    read.reference_start = read.next_reference_start = position
    if not read.is_unmapped:
        read.cigarstring = '20M'
    read.set_tags([('NH', 1), ('HI', 1), ('GX', 'g1'), ('CB', 'pool1'),
                   ('XP', 'pool1'), ('XI', 'original-cell'), ('UB', umi), ('ZI', identity)])
    return read


def test_fresh_bulk_qc_counts_complete_and_singleton_templates_end_to_end(tmp_path, monkeypatch):
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
    (raw / 'matrix.mtx').write_text('%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 3\n')
    (solo / 'CellReads.stats').write_text('featureU\tcountedU\n4\t4\n')
    prep = tmp_path / 'preparation.json'
    prep.write_text(json.dumps(dict(routed_pairs=4, retained_pairs=4,
        barcodes={tag: {'ACGT': 4} for tag in ('L1', 'L2', 'L3', 'SB')})))
    (tmp_path / 's_g.Log.final.out').write_text('Uniquely mapped reads number |\t4\nNumber of reads mapped to multiple loci |\t0\n')
    source = tmp_path / 'aligned.by_name.bam'
    with pysam.AlignmentFile(source, 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as bam:
        # One duplicate input pair, plus both orientations of retained singleton.
        templates = [('a_complete', 'AAAAAAAAAA', None), ('b_duplicate', 'AAAAAAAAAA', None),
                     ('c_r1_only', 'CCCCCCCCCC', 2), ('d_r2_only', 'GGGGGGGGGG', 1)]
        for identity, (name, umi, unmapped_mate) in enumerate(templates, 1):
            for mate in (1, 2):
                flag = 0x1 | (0x40 if mate == 1 else 0x80)
                if mate == unmapped_mate:
                    flag |= 0x4
                elif unmapped_mate:
                    flag |= 0x8
                bam.write(make_record(bam.header, name, flag, 20 * identity, umi, identity))
    summarize('s_g', index, prep, source, solo)
    final_bam = tmp_path / 's_g.rna.dedup.bam'
    pysam.sort('-o', str(final_bam), 'representatives.bam')
    pysam.index(str(final_bam))
    with pysam.AlignmentFile(final_bam) as bam:
        assert bam.has_index()
        records = list(bam.fetch(until_eof=True))
    assert len(records) == 6
    primary_r1 = [r for r in records if r.is_read1 and not r.is_secondary and not r.is_supplementary]
    assert len(primary_r1) == 3
    assert sum(r.is_unmapped for r in primary_r1) == 1
    assert sum(r.is_unmapped for r in records if r.is_read2) == 1
    qc_path = tmp_path / 's_g.rna.qc.json'
    original = json.loads(qc_path.read_text())
    assert original['molecules'] == 3 and original['featureU'] == 4
    spec = dict(restart=True, targets=[dict(sample='s', group='g', modality='rna', id='s_g', library='test')])
    label = 'Final deduplicated BAM read-pair representatives · all barcodes'
    assert f'{label}</span><strong>unavailable' in render(spec, tmp_path)
    digest = hashlib.sha256(final_bam.read_bytes()).hexdigest()
    # Run the same producer CLI invoked by the production module after sorting.
    subprocess.run([sys.executable, str(REPO / 'bin/bulk_rna.py'), '--final-bam-qc',
                    str(final_bam), str(qc_path)], check=True)
    qc = json.loads(qc_path.read_text())
    assert qc['final_bam_pairs'] == len(primary_r1) == 3
    assert qc['final_bam_pairs'] != qc['featureU']
    assert qc['final_bam_unmapped_r1_records'] == 1
    assert qc['final_bam_primary_r2_records'] == 3
    assert qc['final_bam_pairs_match_molecules'] is True
    assert qc['final_bam_pairs_required_flags'] == 0x40
    assert qc['final_bam_pairs_excluded_flags'] == 0x900
    assert qc['final_bam_pairs_source'].startswith(final_bam.name + ': sequential BAM scan')
    assert 'unmapped R1' in qc['final_bam_pairs_definition']
    assert {k: qc[k] for k in original} == original
    assert hashlib.sha256(final_bam.read_bytes()).hexdigest() == digest
    output = tmp_path / 'report'
    write_outputs(spec, tmp_path, output)
    html = (output / 'bulk_report.html').read_text()
    assert f'{label}</span><strong>3</strong>' in html
    assert qc['final_bam_pairs_source'] in html
    assert 'require READ1 (0x40)' in html
    assert 'SUPPLEMENTARY (0x800)' in html
    assert (output / 'rna.gene_by_group.tsv').read_text() == 'gene_id\ts_g\ng1\t3\n'


def test_counting_flags_and_molecule_comparison_are_independent(tmp_path):
    bam_path = tmp_path / 'flags.rna.dedup.bam'
    flags = [0x41, 0x45, 0x40, 0x41 | 0x400, 0x41 | 0x200,
             0x41 | 0x100, 0x41 | 0x800, 0x81]
    with pysam.AlignmentFile(bam_path, 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as bam:
        for i, flag in enumerate(flags):
            bam.write(make_record(bam.header, f'read{i}', flag, i * 20, 'AAAAAAAAAA', i))
    # No index: the audit must count through EOF and include unmapped records.
    qc_path = tmp_path / 'flags.rna.qc.json'
    qc_path.write_text(json.dumps(dict(molecules=999, featureU=1234)))
    update_final_bam_qc(bam_path, qc_path)
    qc = json.loads(qc_path.read_text())
    assert qc['final_bam_pairs'] == 5  # Not 999 molecules or 1234 assigned pairs.
    assert qc['final_bam_unmapped_r1_records'] == 1
    assert qc['final_bam_primary_r2_records'] == 1
    assert qc['final_bam_pairs_match_molecules'] is False


def test_audit_is_wired_after_final_bam_creation_and_before_coverage():
    module = (REPO / 'modules/local/bulk_rna/main.nf').read_text()
    audit = module.index('--final-bam-qc')
    assert module.index('index "${splitName}.rna.dedup.bam"') < audit
    assert audit < module.index('collate -u -O')
    assert '"${splitName}.rna.dedup.bam" "${splitName}.rna.qc.json"' in module[audit:]
    assert 'path("${splitName}.rna.qc.json"), emit: qc' in module
    rna = (REPO / 'subworkflows/local/rna_core/main.nf').read_text()
    assert 'ch_bulk_qc = BULK_RNA_MOLECULES.out.qc' in rna
    workflow = (REPO / 'workflows/treseq.nf').read_text()
    assert 'RNA_CORE.out.bulk_counts.mix(RNA_CORE.out.bulk_qc' in workflow
