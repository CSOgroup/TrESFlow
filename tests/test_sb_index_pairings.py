"""Physical SB routing and explicitly declared logical partition identities."""
import copy
import csv
import json
import os
import shutil
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from test_full_cell_identity import REPO, LOOKUP, LIGATIONS, group, parse_cases, sheet, split_fixture, run, nextflow_test_env
from tresflow_fastq_utils import load_sb_identities


def records(result):
    assert 'error' not in result, result
    return [r for row in result['samples'] for r in row['cell_identity_records']]


def test_available_and_unavailable_chemistries(tmp_path):
    groups = {'B':group(sb_oligo_indices=[13,14,15,16])}
    cases = [sheet(groups, dna=None), sheet(groups, dna='single'),
             sheet({'B':group(sb_oligo_indices=[7])}), sheet(groups),
             sheet({'B':group(sb_oligo_pairings=[{'sb_index':13,'rna_oligo_index':13,'dna_oligo_index':7}])}),
             sheet({'B':group(sb_oligo_pairings=[{'sb_index':7,'rna_oligo_index':7,'dna_oligo_index':13}])})]
    results = parse_cases(tmp_path, cases)
    for result in results[:3]:
        for r in records(result):
            assert r['sb_index'] == r['oligo_index']
    expected = dict(zip(('13','14','15','16'), ('GATC','TGAC','CATG','TACG')))
    for result in results[:2]:
        assert {r['oligo_index']:r['sb_bc'] for r in records(result)} == expected
    for i in (3,5):
        assert 'no available dual-DNA sequence' in results[i]['error']
        assert "physical oligo '13'" in results[i]['error']
    assert {(r['modality'],r['oligo_index'],r['sb_index'],r['sb_bc']) for r in records(results[4])} == {
        ('rna','13','13','GATC'),('dna','07','13','GAC')}
    assert len(results[4]['warnings']) == 1


def test_independent_selections_and_legacy_precedence(tmp_path):
    cases = [
        sheet({'B':group(rna_sb_oligo_indices=[7,13],dna_sb_oligo_indices=[10,7])}),
        sheet({'B':group(rna_sb_barcodes=['GCTA','GATC'],dna_sb_barcodes=['GCA','GAC'])}),
        sheet({'B':group(rna_sb_oligo_indices=[7],dna_sb_barcodes=['GCA'])}),
        sheet({'B':group(rna_sb_barcodes=['GCTA'],dna_sb_oligo_indices=[10])}),
        sheet({'R':{'rna_sb_oligo_indices':[13]},'D':group(dna_sb_oligo_indices=[7])}),
        sheet({'R':{'rna_sb_barcodes':['GCAT']},'D':group(dna_sb_barcodes=['GAT'])}),
        sheet({'B':group(sb_barcodes=['GCAT'],rna_sb_barcodes=['GCTA'],dna_sb_barcodes=['TGCA'])},dna='single'),
        sheet({'B':group(sb_barcodes=['GCAT'],rna_sb_barcodes=['GCTA'])},dna='single'),
        sheet({'B':group(rna_sb_oligo_indices=[7,1],dna_sb_oligo_indices=[1,7])}),
        sheet({'B':group(sb_barcodes=['GCTA'],dna_sb_oligo_indices=[10])}),
    ]
    results = parse_cases(tmp_path, cases)
    for result in results:
        assert all(r['sb_index']==r['oligo_index'] for r in records(result))
    for i in (0,1,2,3,6,7,9):
        assert len(results[i]['warnings'])==1
        assert 'full identifiers will differ' in results[i]['warnings'][0]
        assert 'list order' in results[i]['warnings'][0]
    assert not results[4]['warnings'] and not results[8]['warnings']
    assert 'group namespaces differ' in results[5]['warnings'][0]
    assert {r['modality']:r['oligo_index'] for r in records(results[7])} == {'rna':'07','dna':'01'}


@pytest.mark.parametrize('dna', ['single','dual'])
def test_pairings_resolve_indices_sequences_and_modality_partitions(tmp_path,dna):
    dna_sequence = 'TGCA' if dna=='single' else 'GCA'
    entries = [
        {'sb_index':'007','rna_oligo_index':'07','dna_oligo_index':10},
        {'sb_index':100,'rna_sb_barcode':'GATC'},
        {'sb_index':'100','dna_sb_barcode':'GAC' if dna=='dual' else 'CATG'},
        {'sb_index':'999','rna_oligo_index':'16'},
        {'sb_index':'888','dna_oligo_index':'02'},
    ]
    by_index = sheet({'B':group(sb_oligo_pairings=entries)},dna=dna)
    by_sequence = copy.deepcopy(by_index)
    by_sequence['samples']['VTD11_VTR12']['groups']['B']['sb_oligo_pairings'][0] = {
        'sb_index':'07','rna_sb_barcode':'GCTA','dna_sb_barcode':dna_sequence}
    results = parse_cases(tmp_path,[by_index,by_sequence])
    def resolved(result):
        return {(r['modality'],r['oligo_index'],r['sb_bc'],r['sb_index']) for r in records(result)}
    assert resolved(results[0])==resolved(results[1])
    assert ('dna','10',dna_sequence,'07') in resolved(results[0])
    assert ('rna','13','GATC','100') in resolved(results[0])
    assert len(results[0]['warnings'])==4
    info = tmp_path/'out0/pipeline_info'
    audit = list(csv.DictReader((info/'sb_physical_to_logical.tsv').open(),delimiter='\t'))
    assert {(r['modality'],r['oligo_index'],r['sb_bc'],r['sb_index']) for r in audit}==resolved(results[0])
    warning = (info/'sb_identity_warnings.txt').read_text()
    assert 'SB IDENTITY WARNING' in warning and "Sample 'VTD11_VTR12', group 'B'" in warning
    assert f'actual DNA oligo 10 ({dna_sequence})' in warning
    assert f'DNA chemistry {dna}' in warning and 'shared sb_index 07' in warning
    assert 'one biological partition' in warning
    for row in results[0]['samples']:
        derived = list(csv.DictReader(Path(row['sb_group_map']).open(),delimiter='\t'))
        assert list(derived[0])[-1]=='sb_index'
        assert row['sb_mapping_sha256'] and row['cell_id_version']=='full-cell-v2'


def test_conflicting_selectors_and_ambiguous_routing(tmp_path):
    valid = {'sb_index':7,'rna_oligo_index':7,'dna_oligo_index':10}
    cases=[]
    for field in ('sb_oligo_indices','rna_sb_oligo_indices','dna_sb_oligo_indices',
                  'sb_barcodes','rna_sb_barcodes','dna_sb_barcodes'):
        cases.append(sheet({'B':group(sb_oligo_pairings=[valid],**{field:[7]})}))
    for field in ('rna_sb_oligo_indices','dna_sb_oligo_indices'):
        cases.append(sheet({'B':group(sb_oligo_indices=[7],**{field:[7]})}))
    for modality in ('rna','dna'):
        for sequence_field in ('sb_barcodes',f'{modality}_sb_barcodes'):
            cases.append(sheet({'B':group(**{f'{modality}_sb_oligo_indices':[7],sequence_field:['GCTA']})},
                               dna='single' if modality=='dna' and sequence_field=='sb_barcodes' else 'dual'))
    invalid_indices = [True,False,0,-1,1.0,'+1','1.0',' 1','1 ','','０１',None]
    for index in invalid_indices:
        cases.append(sheet({'B':group(sb_oligo_pairings=[{**valid,'sb_index':index}])}))
        cases.append(sheet({'B':group(sb_oligo_pairings=[{**valid,'dna_oligo_index':index}])}))
        cases.append(sheet({'B':group(rna_sb_oligo_indices=[index],dna_sb_oligo_indices=[10])}))
    for entry in ({}, {'sb_index':7}, {**valid,'cell_index':7},
                  {**valid,'rna_sb_barcode':'GCTA'}, {**valid,'dna_sb_barcode':'GCA'},
                  {**valid,'dna_oligo_index':99}, {'sb_index':7,'rna_sb_barcode':'AAAA'},
                  {'sb_index':7,'dna_sb_barcode':'-'}, {'sb_index':7,'dna_sb_barcode':'TGCA'}):
        cases.append(sheet({'B':group(sb_oligo_pairings=[entry])}))
    for value in ([],None,{},'07',[None]):
        cases.append(sheet({'B':group(sb_oligo_pairings=value)}))
    cases += [
        sheet({'B':group(sb_oligo_pairings=[valid,{**valid,'sb_index':8}])}),
        sheet({'B':group(sb_oligo_pairings=[valid,{'sb_index':'07','rna_oligo_index':8}])}),
        sheet({'B':group(sb_oligo_pairings=[valid]),'C':group(sb_oligo_pairings=[{**valid,'sb_index':8}])}),
        sheet({'B':group(sb_oligo_pairings=[valid])},rna=False),
        sheet({'B':group(sb_oligo_pairings=[valid])},dna=None),
    ]
    results=parse_cases(tmp_path,cases)
    assert all('error' in result for result in results), [(i,r) for i,r in enumerate(results) if 'error' not in r]


@pytest.mark.parametrize('dna', ['single','dual'])
@pytest.mark.parametrize('sequences', [False,True])
def test_real_and_python_split_shared_identity_and_physical_tags(tmp_path,dna,sequences):
    dna_sequences = ('TGCA','CATG') if dna=='single' else ('GCA','GAC')
    dna_indices = ('10','15') if dna=='single' else ('10','07')
    entries = [{'sb_index':'07','rna_oligo_index':7,'dna_oligo_index':10},
               {'sb_index':'13','rna_oligo_index':13,'dna_oligo_index':int(dna_indices[1])}]
    if sequences:
        entries = [{'sb_index':label,'rna_sb_barcode':rna,'dna_sb_barcode':dna_sb}
                   for label,rna,dna_sb in zip(('07','13'),('GCTA','GATC'),dna_sequences)]
    result = parse_cases(tmp_path/'contract',[sheet({'K422_B':group(sb_oligo_pairings=entries)},dna=dna)])[0]
    records(result)
    expected = {f'VTD11_VTR12_K422_B_{i}_{LIGATIONS}' for i in ('07','13')}
    cells={}
    for row in result['samples']:
        modality=row['modality']
        physical = ('GCTA','GATC') if modality=='rna' else dna_sequences
        indices = ('07','13') if modality=='rna' else dna_indices
        outputs={mode:split_fixture(tmp_path/modality/mode,modality,mode,physical,indices,
                                   sb_map=row['sb_group_map'],mo_map=row.get('mo_map')) for mode in ('mock','real')}
        stem='VTD11_VTR12_K422_B'+('_H3K27ac' if modality=='dna' else '')
        for mate in ('R1','R2'):
            text=(outputs['real']/f'{stem}_{mate}.fastq').read_text()
            assert text==(outputs['mock']/f'{stem}_{mate}.fastq').read_text()
            headers=text.splitlines()[::4]
            tags=[dict(t.split(':',2)[::2] for t in h.split()[1:]) for h in headers]
            cells[modality]={t['CB'] for t in tags}
            assert cells[modality]==expected
            assert all(t['CB']==t['XI'] and t['UM']=='ACGTACGTAC' for t in tags)
            assert {t['SB'] for t in tags}==set(physical)
            assert {t['MO'] for t in tags}=={'GCCTCTAT'}
            if modality=='dna': assert {t['RG'] for t in tags}=={'I:R:F:L1'}
        if modality=='rna':
            for mode in ('mock','real'):
                out=outputs[mode]
                run([sys.executable,REPO/'bin/run_fq_to_sam.py','--mode',mode,
                     '--script',REPO/'scripts/core_runtime/FqToSAM.codon','--r1',out/f'{stem}_R1.fastq',
                     '--r2',out/f'{stem}_R2.fastq','--output-sam',out/'rna.sam'])
            assert (outputs['real']/'rna.sam').read_bytes()==(outputs['mock']/'rna.sam').read_bytes()
            solo=tmp_path/'Solo.outGeneFull'
            run([sys.executable,REPO/'bin/mock_starsolo_outputs.py',outputs['real']/'rna.sam',solo])
            assert set((solo/'filtered/barcodes.tsv').read_text().splitlines())==expected
    assert cells['rna']==cells['dna']


@pytest.mark.parametrize('modality',['rna','dna'])
def test_new_indices_python_and_codon_split(tmp_path,modality):
    result=parse_cases(tmp_path/'contract',[sheet({'K422_B':group(sb_oligo_indices=[13,14,15,16])},dna='single')])[0]
    row=next(r for r in result['samples'] if r['modality']==modality)
    outputs={mode:split_fixture(tmp_path/mode,modality,mode,('GATC','TGAC','CATG','TACG'),('13','14','15','16'),
                               sb_map=row['sb_group_map'],mo_map=row.get('mo_map')) for mode in ('mock','real')}
    for path in outputs['real'].glob('*.fastq'):
        assert path.read_bytes()==(outputs['mock']/path.name).read_bytes()
        assert all(f'CB:Z:VTD11_VTR12_K422_B_{i}_' in path.read_text() for i in ('13','14','15','16'))


def test_warnings_are_logged_before_preflight_and_invalid_sheets_do_not_warn(tmp_path):
    assert shutil.which('nextflow')
    cases=[sheet({'B':group(sb_oligo_pairings=[{'sb_index':7,'rna_oligo_index':7,'dna_oligo_index':10}])}),
           sheet({'B':group(rna_sb_oligo_indices=[7],dna_sb_oligo_indices=[10])})]
    invalid=copy.deepcopy(cases[0]); invalid['samples']['VTD11_VTR12']['groups']['B']['sb_oligo_indices']=[7]
    for i,case in enumerate(cases+[invalid]):
        root=tmp_path/str(i);root.mkdir()
        # Samplesheet parsing succeeds, then this intentionally absent runtime
        # environment proves warnings were emitted before any preflight work.
        case['runtime']['env_prefix']=str(root/'absent_environment')
        path=root/'sheet.yaml';path.write_text(json.dumps(case))
        result=subprocess.run(['nextflow','run',str(REPO),'-preview','-ansi-log','false',
                               '--samplesheet',str(path),'--outdir',str(root/'out')],cwd=root,
                              env=nextflow_test_env(root/'nxf_home'),
                              capture_output=True,text=True,timeout=90)
        assert result.returncode!=0
        log=(root/'.nextflow.log').read_text()
        if i<2:
            assert log.index('SB IDENTITY WARNING') < log.index('TrESFlow storage paths:')
            assert 'SB IDENTITY WARNING' in result.stdout+result.stderr
            assert (root/'out/pipeline_info/sb_physical_to_logical.tsv').is_file()
            assert 'SB IDENTITY WARNING' in (root/'out/pipeline_info/sb_identity_warnings.txt').read_text()
        else:
            assert 'SB IDENTITY WARNING' not in log
            assert not (root/'out/pipeline_info/sb_identity_warnings.txt').exists()


@pytest.mark.parametrize('paired',[False,True])
def test_bounded_workflow_accepts_mismatches_and_preserves_identifiers(tmp_path,paired):
    import yaml
    groups={}
    for name,rna,dna in [('Normal','05','01'),('Co2','06','02')]:
        selector={'sb_oligo_pairings':[{'sb_index':rna,'rna_oligo_index':rna,'dna_oligo_index':dna}]} if paired else {
            'rna_sb_oligo_indices':[rna],'dna_sb_oligo_indices':[dna]}
        groups[name]={'mark_barcodes':{'M':'AGGCTATA'},**selector}
    case=sheet(groups,sample='workflow_sample')
    case['runtime']={'env_prefix':os.environ.get('CONDA_PREFIX',str(Path(shutil.which('STAR')).parent.parent)),
                     'tmpdir':str(tmp_path/'runtime_tmp')}
    case['samples']['workflow_sample']['dna']['reads']['i1']=str(REPO/'assets/testdata/test_dna_I1_dual_lig.fastq')
    path=tmp_path/'sheet.yaml';path.write_text(yaml.safe_dump(case,sort_keys=False))
    result=subprocess.run(['nextflow','run',str(REPO),'-profile','test','-ansi-log','false',
                           '-work-dir',str(tmp_path/'work'),'--samplesheet',str(path),'--outdir',str(tmp_path/'out'),
                           '--cleanup_work','false','--publish_split_fastqs','true'],cwd=tmp_path,
                          env=nextflow_test_env(tmp_path/'nxf_home'),
                          capture_output=True,text=True,timeout=180)
    (tmp_path/'run.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stdout+result.stderr
    info=tmp_path/'out/pipeline_info'
    warning=(info/'sb_identity_warnings.txt').read_text()
    log=(tmp_path/'.nextflow.log').read_text()
    assert warning.count('SB IDENTITY WARNING')==2
    assert log.index('SB IDENTITY WARNING')<log.index('TrESFlow storage paths:')
    audit=list(csv.DictReader((info/'sb_physical_to_logical.tsv').open(),delimiter='\t'))
    assert {r['sb_index'] for r in audit if r['modality']=='dna'}==({'05','06'} if paired else {'01','02'})
    for name,rna,dna,rna_sb,dna_sb in [('Normal','05','01','CAGT','GAT'),('Co2','06','02','ACGT','AGT')]:
        for modality,physical,sb in [('rna',rna,rna_sb),('dna',dna,dna_sb)]:
            logical=rna if paired else physical
            stem=f'workflow_sample_{name}'+('_M' if modality=='dna' else '')
            with gzip.open(tmp_path/f'out/{modality}_split_fastqs/{stem}_R1.fastq.gz','rt') as handle:
                text=handle.read()
            if modality=='dna' and name=='Co2':
                assert text==''  # existing raw DNA fixture contains only physical 01
                continue
            assert f'CB:Z:workflow_sample_{name}_{logical}_' in text
            assert f'XI:Z:workflow_sample_{name}_{logical}_' in text
            assert f'SB:Z:{sb}' in text
            if modality=='rna':
                solo=tmp_path/f'out/rna_align/{stem}.Solo.outGeneFull'
                assert f'workflow_sample_{name}_{logical}_' in (solo/'raw/barcodes.tsv').read_text()


def test_schema_and_existing_examples(tmp_path):
    import jsonschema
    import yaml
    schema=json.loads((REPO/'assets/samplesheet.schema.json').read_text())
    jsonschema.Draft202012Validator.check_schema(schema)
    paths=[REPO/'assets/samplesheet.example.yaml',REPO/'assets/samplesheet.sb-pairings.example.yaml',
           REPO/'tests/samplesheets/group_specific_modalities.yaml',
           REPO/'tests/samplesheets/multi_dna_single_sequence.yaml',
           REPO/'tests/samplesheets/multi_dna_dual_sequence.yaml',
           REPO/'tests/samplesheets/multi_rna_comma.yaml']
    cases=[]
    for path in paths:
        case=yaml.safe_load(path.read_text())
        jsonschema.validate(case,schema)
        # Rebase relative file paths when copying the legacy sheets to temp.
        def absolute(value):
            if isinstance(value,list): return [absolute(v) for v in value]
            if isinstance(value,str):
                return str((path.parent/value).resolve())
            return value
        case['references']={k:absolute(v) if k!='species' and isinstance(v,str) else v for k,v in case['references'].items()}
        for sample in case['samples'].values():
            for modality in ('rna','dna'):
                if modality in sample:
                    for role,value in sample[modality]['reads'].items():
                        if isinstance(value,str) and ',' in value: # established comma-list interface
                            value=','.join(absolute(v.strip()) for v in value.split(','))
                        else: value=absolute(value)
                        sample[modality]['reads'][role]=value
        cases.append(case)
    results=parse_cases(tmp_path,cases)
    assert all('error' not in r for r in results),results
    for i,result in enumerate(results):
        if i!=1: assert all(r['oligo_index']==r['sb_index'] for r in records(result))
    invalid=sheet({'B':group(sb_oligo_pairings=[{'sb_index':True,'rna_oligo_index':7}])})
    with pytest.raises(jsonschema.ValidationError): jsonschema.validate(invalid,schema)


def test_derived_logical_selector_is_explicit_and_strict(tmp_path):
    path=tmp_path/'map.tsv'
    path.write_text('sample\tsb_group\tsb_bc\toligo_index\tsb_index\nS\tG\tGCTA\t07\t13\n')
    record=load_sb_identities(path,'S')['GCTA']
    assert (record.oligo_index,record.sb_bc,record.sb_index)==('07','GCTA','13')
    path.write_text('sample\tsb_group\tsb_bc\toligo_index\nS\tG\tGCTA\t07\n')
    assert load_sb_identities(path,'S')['GCTA'].sb_index=='07'
    for label in ('', '0','-1','+1','1.0',' 1','1 ','１２'):
        path.write_text(f'sample\tsb_group\tsb_bc\toligo_index\tsb_index\nS\tG\tGCTA\t07\t{label}\n')
        with pytest.raises(ValueError,match='sb_index'): load_sb_identities(path,'S')


@pytest.mark.parametrize('label', [' 13', '13 ', ''])
def test_codon_rejects_malformed_derived_logical_indices(tmp_path, label):
    path = tmp_path/'map.tsv'
    path.write_text(f'sample\tsb_group\tsb_bc\toligo_index\tsb_index\nS\tG\tGCTA\t07\t{label}\n')
    with pytest.raises(ValueError, match='sb_index'):
        load_sb_identities(path, 'S')
    for mate in ('r1', 'r2'):
        (tmp_path/f'{mate}.fastq').write_text(
            f'@I:R:F:1:1:1:1:U CB:Z:GCTA{LIGATIONS}\tSB:Z:GCTA\tUM:Z:ACGTACGTAC\n'
            'ACGTACGT\n+\nIIIIIIII\n')
    out = tmp_path/'out'; out.mkdir()
    result = subprocess.run([
        'codon', 'run', '-plugin', 'seq', '-release',
        str(REPO/'scripts/core_runtime/Split_ReadsV2.codon'), 'S', str(out), 'TEST', 'rna', '-',
        str(tmp_path/'r1.fastq'), str(tmp_path/'r2.fastq'), str(path)],
        capture_output=True, text=True, timeout=90)
    assert result.returncode != 0
    assert 'Missing or invalid sb_index in derived SB map' in result.stderr
    assert not list(out.glob('*.fastq'))


def test_schema_preserves_legacy_chemistry_and_sequence_normalization(tmp_path):
    import jsonschema
    schema = json.loads((REPO/'assets/samplesheet.schema.json').read_text())
    cases = []
    for chemistry in ('single', 'dual'):
        for indexed in (True, False):
            selectors = {'sb_oligo_indices':[7]} if indexed else {
                'sb_barcodes':[' GCTA '],
                'dna_sb_barcodes':[' GAC ' if chemistry == 'dual' else ' GCTA ']}
            for spelling in (chemistry, chemistry.upper(), f'\t{chemistry.title()} \n'):
                case = sheet({'B':group(**selectors)}, dna=chemistry)
                case['samples']['VTD11_VTR12']['dna']['tagmentation'] = spelling
                jsonschema.validate(case, schema)
                cases.append(case)
    results = parse_cases(tmp_path, cases)
    for offset in range(0, len(results), 3):
        assert records(results[offset]) == records(results[offset+1]) == records(results[offset+2])
        assert all(r['oligo_index'] == r['sb_index'] == '07' for r in records(results[offset]))
    # Chemistry normalization must also apply to the single-DNA conflict rule.
    conflict = sheet({'B':group(sb_barcodes=['GCTA'], dna_sb_oligo_indices=[7])}, dna='single')
    conflict['samples']['VTD11_VTR12']['dna']['tagmentation'] = ' Single '
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(conflict, schema)
    assert 'conflicting selectors' in parse_cases(tmp_path/'conflict', [conflict])[0]['error']
