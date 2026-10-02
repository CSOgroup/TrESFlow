"""Exercise actual Nextflow hashes and publication with a bounded mock graph."""
import csv
import json
import os
import subprocess
from pathlib import Path

import pytest

from test_full_cell_identity import nextflow_test_env

REPO = Path(__file__).resolve().parents[1]


def test_lookup_content_and_identity_changes_invalidate_resumed_tasks(tmp_path):
    import shutil
    import yaml
    assert shutil.which('nextflow'), 'Real Nextflow required'
    data = REPO / 'assets/testdata'
    root = data / 'TrESFlow_References'
    sheet = {
        'library_name': 'CACHE_TEST',
        'runtime': {'env_prefix': os.environ.get('CONDA_PREFIX', str(Path(shutil.which('STAR')).parent.parent)),
                    'tmpdir': str(tmp_path / 'runtime_tmp')},
        'references': {'species':'human', 'root':str(root), 'ligation_barcode_whitelist':str(root/'ligation_barcode_whitelist.txt'),
                       'rna_ref_dir':str(root/'rna/human/star'), 'dna_ref_dir':str(root/'dna/human/bwa'),
                       'dna_blacklist_bed':str(root/'dna/human/blacklist.bed'), 'dna_chrom_sizes':str(root/'dna/human/chrom.sizes'),
                       'dna_effective_genome_size':12},
        'samples': {'cache_sample': {
            'groups': {'Normal': {'rna_sb_barcodes':['CAGT'], 'dna_sb_barcodes':['GTC'], 'mark_barcodes':{'M':'AGGCTATA'}},
                       'K422_A': {'rna_sb_barcodes':['GCAT'], 'dna_sb_barcodes':['GAT'], 'mark_barcodes':{'M':'AGGCTATA'}}},
            'rna':{'reads':{r:str(data/f'test_rna_{r.upper()}.fastq') for r in ('i1','r1','r2')}},
            'dna':{'tagmentation':'dual','reads':{'i1':str(data/'test_dna_I1_dual_lig.fastq'),
                    'r1':str(data/'test_dna_R1.fastq'), 'r2':str(data/'test_dna_R2.fastq')}}}}
    }
    samplesheet=tmp_path/'samples.yaml';samplesheet.write_text(yaml.safe_dump(sheet,sort_keys=False))
    lookup=tmp_path/'lookup.tsv';lookup.write_bytes((REPO/'assets/sb_oligo_lookup.v1.tsv').read_bytes())
    command=['nextflow','run',str(REPO),'-profile','test','-ansi-log','false','-work-dir',str(tmp_path/'work'),
             '--samplesheet',str(samplesheet),'--outdir',str(tmp_path/'out'), '--sb_oligo_lookup',str(lookup),
             '--cleanup_work','false','--publish_split_fastqs','true']
    env=nextflow_test_env(tmp_path/'nxf_home')
    def execute(number,resume=False):
        trace=tmp_path/f'trace{number}.tsv'
        result=subprocess.run(command+['-with-trace',str(trace)]+(['-resume'] if resume else []),
                              cwd=tmp_path,env=env,capture_output=True,text=True,timeout=180)
        (tmp_path/f'run{number}.log').write_text(result.stdout+result.stderr)
        assert result.returncode==0,result.stdout+result.stderr
        return list(csv.DictReader(trace.open(),delimiter='\t'))
    first=execute(1)
    second=execute(2,True)
    def rows(stage,trace):return [r for r in trace if r['name'].split(' (')[0].endswith(':'+stage)]
    stages=['TAG_RNA_SAMPLE_BARCODE','TAG_DNA_SAMPLE_BARCODE','SPLIT_RNA_READS','SPLIT_DNA_READS',
            'FQ_TO_SAM','RNA_STARSOLO_ALIGN','ALIGN_DNA','GATK4_MARKDUPLICATES']
    downstream_stages = [
        'TAG_RNA_UMI', 'TAG_RNA_CELL_BARCODE', 'TAG_DNA_MODALITY_BARCODE', 'TAG_DNA_CELL_BARCODE',
        'TRIM_RNA_FASTQS', 'TRIM_DNA_FASTQS', 'DUAL_TAG_ARTIFACT_FILTER', 'BARCODE_GATE_METRICS',
        'COMPRESS_RNA_SPLIT_FASTQS', 'COMPRESS_DNA_SPLIT_FASTQS', 'RNA_FILTERED_BAM', 'RNA_COVERAGE',
        'FILTER_CANONICAL_DNA_ALIGNED_BAM', 'NORMALIZE_DNA_MARKDUPLICATES', 'SPLIT_DUPLICATES_DNA',
        'DEEPTOOLS_BAMCOVERAGE', 'TRES_REPORT_HTML', 'MULTIQC',
    ]
    for stage in stages:
        assert rows(stage,first) and all(r['status']=='CACHED' for r in rows(stage,second)),(stage,second)
    # Change content at the same path, preserving byte size and mtime. The
    # explicit SHA-256 in metadata must invalidate even standard path caching.
    stat=lookup.stat()
    lookup.write_text(lookup.read_text().replace('05\tCAGT\tGTC','17\tCAGT\tGTC'))
    os.utime(lookup,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    assert lookup.stat().st_size==stat.st_size
    third=execute(3,True)
    for stage in stages + downstream_stages:
        assert rows(stage,first) and rows(stage,third), stage
        assert all(r['status']=='COMPLETED' for r in rows(stage,third)),(stage,third)
        assert {r['hash'] for r in rows(stage,first)}.isdisjoint({r['hash'] for r in rows(stage,third)})
    import gzip
    fastq=tmp_path/'out/rna_split_fastqs/cache_sample_Normal_R1.fastq.gz'
    with gzip.open(fastq,'rt') as handle:text=handle.read()
    assert 'CB:Z:cache_sample_Normal_17_' in text
    barcodes=(tmp_path/'out/rna_align/cache_sample_Normal.Solo.outGeneFull/raw/barcodes.tsv').read_text()
    assert 'cache_sample_Normal_17_' in barcodes and '_05_' not in barcodes

    # Switch to explicit pairing with unchanged physical sequences, then change
    # only the logical sb_index at the same samplesheet path on another resume.
    normal = sheet['samples']['cache_sample']['groups']['Normal']
    normal.pop('rna_sb_barcodes'); normal.pop('dna_sb_barcodes')
    normal['sb_oligo_pairings'] = [{'sb_index':'42', 'rna_sb_barcode':'CAGT', 'dna_sb_barcode':'GTC'}]
    samplesheet.write_text(yaml.safe_dump(sheet, sort_keys=False))
    fourth = execute(4, True)
    physical_before = list(csv.DictReader((tmp_path/'out/pipeline_info/sb_physical_to_logical.tsv').open(), delimiter='\t'))
    normal['sb_oligo_pairings'][0]['sb_index'] = '43'
    samplesheet.write_text(yaml.safe_dump(sheet, sort_keys=False))
    fifth = execute(5, True)
    for stage in stages + downstream_stages:
        assert rows(stage, fourth) and rows(stage, fifth), stage
        assert all(r['status']=='COMPLETED' for r in rows(stage, fifth)), (stage, fifth)
        assert {r['hash'] for r in rows(stage, fourth)}.isdisjoint({r['hash'] for r in rows(stage, fifth)})
    physical_after = list(csv.DictReader((tmp_path/'out/pipeline_info/sb_physical_to_logical.tsv').open(), delimiter='\t'))
    assert [{k:v for k,v in r.items() if k!='sb_index'} for r in physical_before] == [
        {k:v for k,v in r.items() if k!='sb_index'} for r in physical_after]
    assert {r['sb_index'] for r in physical_after if r['sb_group']=='Normal'} == {'43'}
    for number, label in [(4,'42'), (5,'43')]:
        log = (tmp_path/f'run{number}.log').read_text()
        assert 'SB IDENTITY WARNING' in log and f'shared sb_index {label}' in log
    assert 'shared sb_index 43' in (tmp_path/'.nextflow.log').read_text()
    assert 'shared sb_index 43' in (tmp_path/'out/pipeline_info/sb_identity_warnings.txt').read_text()
    with gzip.open(fastq, 'rt') as handle: text = handle.read()
    assert 'CB:Z:cache_sample_Normal_43_' in text and 'SB:Z:CAGT' in text
    assert 'cache_sample_Normal_43_' in (tmp_path/'out/rna_align/cache_sample_Normal.Solo.outGeneFull/raw/barcodes.tsv').read_text()
