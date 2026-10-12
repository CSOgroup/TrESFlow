"""Behavioral RNA export tests; pigz and FASTQ validation execute for real."""
import csv
import gzip
import importlib.util
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
spec = importlib.util.spec_from_file_location('rna_export', REPO / 'bin/export_rna_fastqs.py')
EXPORT = importlib.util.module_from_spec(spec)
spec.loader.exec_module(EXPORT)
from tresflow_fastq_utils import fastq_iter, parse_header

UMI = 'GCTACCGGTT'  # Actual reverse-complemented tag; not raw AACCGGTAGC.
CB = 'sample_with_underscores_Pool_A_05_ACGTACGTTGCATGCAGATCGATC'
COMMENT = f'UM:Z:{UMI}\tCB:Z:{CB}\tXI:Z:{CB}\tSB:Z:CAGT\tRG:Z:lane1'


def record(name, sequence, comment=COMMENT):
    return (f'@{name} {comment}', sequence, '+', ''.join(chr(33 + i % 40) for i in range(len(sequence))))


def test_clip_once_and_copy_oriented_umi_to_identical_names():
    r1, r2 = record('original_id/1', 'G' * 40), record('original_id/2', 'ACTGAACCGGTAGC' + 'ACGT' * 10)
    pair, pool, tags = EXPORT.convert_pair(r1, r2, 14, 10, 20)
    assert pair is not None
    one, two = [text.splitlines() for text in pair]
    name = f'original_id_{pool}_{UMI}'
    assert one[0] == two[0] == f'@{name} {COMMENT}'
    assert '_' not in pool and pool.isalnum()
    assert one[1:] == [r1[1], '+', r1[3]]
    assert two[1:] == [r2[1][14:], '+', r2[3][14:]]
    assert two[1] == 'ACGT' * 10  # No second 14-base removal.
    assert tags['UM'] == UMI
    assert UMI != 'AACCGGTAGC'


@pytest.mark.parametrize('r1_length,r2_length,kept', [(20, 34, True), (19, 40, False), (40, 33, False), (40, 14, False)])
def test_minimum_length_applies_to_both_mates_after_clipping(r1_length, r2_length, kept):
    pair, _, _ = EXPORT.convert_pair(record('read', 'A' * r1_length), record('read', 'C' * r2_length), 14, 10, 20)
    assert (pair is not None) == kept


def test_pool_ids_keep_distinct_barcode_provenance_and_are_stable():
    first = EXPORT.pool_identifier(CB, 'CAGT')
    assert first == EXPORT.pool_identifier(CB, 'CAGT')
    assert first != EXPORT.pool_identifier(CB.replace('Pool_A', 'Pool_B'), 'CAGT')
    assert first != EXPORT.pool_identifier(CB, 'GTCA')


@pytest.mark.parametrize('mutation', ['name', 'umi', 'cb', 'missing'])
def test_invalid_pair_identity_fails(mutation):
    r1, r2 = record('read', 'A' * 40), record('other' if mutation == 'name' else 'read', 'C' * 40)
    if mutation == 'umi':
        r2 = record('read', 'C' * 40, COMMENT.replace(UMI, 'AAAAAAAAAA'))
    elif mutation == 'cb':
        r2 = record('read', 'C' * 40, COMMENT.replace('XI:Z:' + CB, 'XI:Z:wrong'))
    elif mutation == 'missing':
        r2 = record('read', 'C' * 40, COMMENT.replace('UM:Z:' + UMI, ''))
    with pytest.raises(ValueError):
        EXPORT.convert_pair(r1, r2, 14, 10, 20)


def write_records(path, records):
    path.write_text(''.join('\n'.join(r) + '\n' for r in records))


def test_gzip_exports_pair_loss_accounting_and_pool_map(tmp_path):
    pigz = shutil.which('pigz')
    if not pigz:
        pytest.skip('pigz unavailable')
    r1, r2 = tmp_path / 'sample_Pool_A_R1.fastq', tmp_path / 'sample_Pool_A_R2.fastq'
    write_records(r1, [record('keep', 'A' * 40), record('short_r2', 'G' * 40), record('short_r1', 'C' * 19)])
    write_records(r2, [record('keep', 'C' * 40), record('short_r2', 'G' * 33), record('short_r1', 'T' * 40)])
    original = (r1.read_bytes(), r2.read_bytes())
    out = tmp_path / 'export'
    EXPORT.export_fastqs([dict(split_name='sample_Pool_A', group='Pool_A', r1=str(r1), r2=str(r2))], out, 'sample', 14, 10, 20, pigz, 1)
    one, two = [list(fastq_iter(out / (p.name + '.gz'))) for p in (r1, r2)]
    assert len(one) == len(two) == 1
    assert parse_header(one[0][0])[0] == parse_header(two[0][0])[0]
    assert one[0][1] == 'A' * 40 and two[0][1] == 'C' * 26
    counts = {r['metric']: int(r['pairs']) for r in csv.DictReader((out / 'sample.rna_export_retention.tsv').open(), delimiter='\t')}
    assert counts == {'export_input_pairs': 3, 'exported_pairs': 1, 'export_short_pairs': 2}
    pools = list(csv.DictReader((out / 'sample.rna_export_pools.tsv').open(), delimiter='\t'))
    assert pools[0]['CB'] == pools[0]['XI'] == CB
    assert pools[0]['SB'] == 'CAGT' and pools[0]['pool_id'] in one[0][0]
    assert (r1.read_bytes(), r2.read_bytes()) == original
    assert not list(out.glob('*.fastq'))


def test_empty_split_exports_valid_empty_gzip_files(tmp_path):
    pigz = shutil.which('pigz')
    if not pigz:
        pytest.skip('pigz unavailable')
    r1, r2 = tmp_path / 's_empty_R1.fastq', tmp_path / 's_empty_R2.fastq'
    r1.touch(); r2.touch()
    out = tmp_path / 'export'
    EXPORT.export_fastqs([dict(split_name='s_empty', group='empty', r1=str(r1), r2=str(r2))], out, 's', 14, 10, 20, pigz, 1)
    assert all(gzip.decompress(p.read_bytes()) == b'' for p in out.glob('*.gz'))


def test_compression_failure_is_not_success(tmp_path):
    r1, r2 = tmp_path / 's_A_R1.fastq', tmp_path / 's_A_R2.fastq'
    r1.touch(); r2.touch()
    with pytest.raises(RuntimeError, match='pigz failed'):
        EXPORT.export_fastqs([dict(split_name='s_A', group='A', r1=str(r1), r2=str(r2))], tmp_path / 'out', 's', 14, 10, 20, '/bin/false', 3)
    assert not (tmp_path / 'out/s.rna_export_retention.tsv').exists()


def test_report_includes_export_losses_without_changing_full_mode(tmp_path):
    sys.path.insert(0, str(REPO / 'lib'))
    from tresflow_qc.core import build_report_model
    from tresflow_qc.html_report import render_report
    (tmp_path / 's.rna_sample_barcode.stats.tsv').write_text('reads\t3\t100\nbc_reads\t3\t100\nreads_without_bc\t0\t0\n')
    fields = 'sample_id\tmodality\tgroup\tbranch\tmetric\tpairs\tunit\n'
    (tmp_path / 's.rna_read_retention.tsv').write_text(fields + 's\trna\t__all__\t__all__\tsplit_input_pairs\t3\tread_pairs\ns\trna\t__all__\t__all__\tjoint_barcode_accepted_pairs\t3\tread_pairs\ns\trna\tA\tRNA\trouted_branch_pairs\t3\tread_pairs\n')
    export = tmp_path / 's.rna_export_retention.tsv'
    export.write_text(fields + 's\trna\tA\tRNA\texport_input_pairs\t3\tread_pairs\ns\trna\tA\tRNA\texported_pairs\t1\tread_pairs\ns\trna\tA\tRNA\texport_short_pairs\t2\tread_pairs\n')
    model = build_report_model(tmp_path, run_metadata={'preprocess_only': True})
    assert model.retention[-1]['stage'] == 'Bulk FASTQ export'
    assert model.retention[-1]['output_pairs'] == 1
    assert next(r['value'] for r in model.qc_metrics if r['metric'] == 'export_short_pairs') == 2
    html = tmp_path / 'report.html'
    render_report(model, [], html)
    assert '3 / 1 (33.3% of raw)' in html.read_text()
    full = build_report_model(tmp_path)
    assert not any(r['stage'] == 'Bulk FASTQ export' for r in full.retention)
    assert not any(r['metric'].startswith('export_') for r in full.qc_metrics)
    export.write_text(export.read_text().replace('export_short_pairs\t2', 'export_short_pairs\t1'))
    with pytest.raises(ValueError, match='routing/export retention mismatch'):
        build_report_model(tmp_path, run_metadata={'preprocess_only': True})


@pytest.mark.parametrize('cpus', [1, 2, 3, 4, 8, 9, 16, 32])
def test_export_cpu_budget_includes_both_mates_and_conversion(cpus):
    threads = EXPORT.compressor_threads(cpus)
    if cpus < 3:
        assert threads == 0
    else:
        actual_threads = 1 if threads == 1 else threads + 2
        assert 1 + 2 * actual_threads <= cpus
    if cpus == 16:
        assert threads == 5


@pytest.mark.parametrize('comment', [
    COMMENT, COMMENT + '\tSB:Z:wrong',
    COMMENT.replace('RG:Z:lane1', 'RG:Z:run:flowcell:lane'),
    COMMENT.replace('UM:Z:' + UMI, 'UM:bad UM:Z:' + UMI),
])
def test_single_pass_tag_parser_preserves_first_match_semantics(comment):
    from tresflow_fastq_utils import find_tag_value
    tags = EXPORT.parse_tags(comment)
    assert {key: tags.get(key, '') for key in EXPORT.REQUIRED_TAGS} == {key: find_tag_value(comment, key) for key in EXPORT.REQUIRED_TAGS}


def test_buffered_cpu_budgets_export_identical_bytes_and_metrics(tmp_path, monkeypatch):
    pigz = shutil.which('pigz')
    if not pigz:
        pytest.skip('pigz unavailable')
    # Small buffer forces multiple flushes and a final partial flush.
    monkeypatch.setattr(EXPORT, 'WRITE_BUFFER_SIZE', 1024)
    r1, r2 = tmp_path / 's_A_R1.fastq', tmp_path / 's_A_R2.fastq'
    rows1 = [record(f'read_{i}', 'A' * 60) for i in range(41)]
    rows2 = [record(f'read_{i}', 'C' * (33 if i % 5 == 0 else 60)) for i in range(41)]
    write_records(r1, rows1); write_records(r2, rows2)
    expected = [[], []]
    for one, two in zip(rows1, rows2):
        pair, _, _ = EXPORT.convert_pair(one, two, 14, 10, 20)
        if pair:
            for destination, text in zip(expected, pair):
                destination.append(text)
    outputs = []
    for cpus in (1, 2, 4, 16):
        out = tmp_path / str(cpus)
        EXPORT.export_fastqs([dict(split_name='s_A', group='A', r1=str(r1), r2=str(r2))], out, 's', 14, 10, 20, pigz, cpus)
        assert [gzip.decompress((out / (p.name + '.gz')).read_bytes()).decode() for p in (r1, r2)] == [''.join(records) for records in expected]
        outputs.append({p.name: p.read_bytes() for p in out.glob('*.tsv')})
    assert all(result == outputs[0] for result in outputs)
