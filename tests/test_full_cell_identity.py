"""Behavioral full-cell-v1 regressions, including real STAR, Codon and GATK.

No private data or external reference downloads are required. All fixture data
are deterministic and generated under pytest's temporary directory.
"""
import csv
import importlib.util
import json
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from test_multi_fastq_inputs import base_contract, find_nextflow_jar, GROOVY_PARSE

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'bin'))
from tresflow_fastq_utils import canonicalize_fastq_comment, canonicalize_dna_fastq_comment

LOOKUP = REPO / 'assets/sb_oligo_lookup.v1.tsv'
LIGATIONS = 'ACGTACGTTGCATGCAGATCGATC'


def run(command, **kwargs):
    result = subprocess.run(list(map(str, command)), capture_output=True, text=True, **kwargs)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def parse_cases(root, cases, lookup=LOOKUP):
    """Validate many independent sheets in one real Groovy JVM."""
    root.mkdir(parents=True, exist_ok=True)
    manifest = []
    for i, sheet in enumerate(cases):
        path = root / f'{i}.yaml'
        path.write_text(json.dumps(sheet))
        manifest.append([str(path), str(root / f'out{i}')])
    (root / 'cases.json').write_text(json.dumps(manifest))
    defaults = GROOVY_PARSE[GROOVY_PARSE.index('def defaults = '):GROOVY_PARSE.index('def contract = ')]
    script = '''
      def parser = new GroovyClassLoader().parseClass(new File(args[0]))
    ''' + defaults + '''
      def cases = new groovy.json.JsonSlurper().parse(new File(args[1]))
      def results = cases.collect { entry ->
        try {
          def contract = parser.parseContract(entry[0], [outdir: entry[1], barcode_defaults: defaults, sb_oligo_lookup: args[2]])
          return [samples: contract.samples]
        } catch(IllegalArgumentException e) { return [error: e.message] }
      }
      println groovy.json.JsonOutput.toJson(results)
    '''
    assert shutil.which('java') and find_nextflow_jar(), 'Real Groovy/Nextflow required'
    output = run(['java', '-cp', find_nextflow_jar(), 'groovy.ui.GroovyMain', '-e', script, '--',
                  REPO / 'lib/SamplesheetParser.groovy', root / 'cases.json', lookup], cwd=REPO)
    return json.loads(output)


def sheet(groups, rna=True, dna='dual', sample='VTD11_VTR12'):
    assets = REPO / 'assets/testdata'
    row = {'groups': groups}
    if rna:
        row['rna'] = {'reads': {role: str(assets / f'test_rna_{role.upper()}.fastq') for role in ('i1', 'r1', 'r2')}}
    if dna:
        row['dna'] = {'tagmentation': dna, 'reads': {role: str(assets / f'test_dna_{role.upper()}.fastq') for role in ('i1', 'r1', 'r2')}}
        if dna == 'single':
            row['dna']['reads']['i2'] = str(assets / 'test_dna_I2_single_lig.fastq')
    return base_contract({sample: row})


def group(**kwargs):
    return {'mark_barcodes': {'H3K27ac': 'GCCTCTAT'}, **kwargs}


def test_index_and_legacy_sequences_resolve_equally_and_keep_underscores(tmp_path):
    cases = [
        sheet({'B': group(sb_oligo_indices=[1, '07']), 'K422_B': group(sb_oligo_indices=['12'])}),
        sheet({'B': group(rna_sb_barcodes=['GCTA', 'GCAT'], dna_sb_barcodes=['GAT', 'GAC']),
               'K422_B': group(rna_sb_barcodes=['TCGA'], dna_sb_barcodes=['TAC'])}),
        sheet({'B': group(sb_oligo_indices=['07'])}, dna='single'),
        sheet({'B': group(sb_barcodes=['GCTA'])}, dna='single'),
        sheet({'A': {'sb_oligo_indices':[1]}}, dna=None),
        sheet({'A': group(sb_oligo_indices=[1])}, rna=False),
        sheet({'RNA_only': {'rna_sb_barcodes':['GACT']},
               'DNA_only': group(dna_sb_barcodes=['ATG'])}),
    ]
    results = parse_cases(tmp_path, cases)
    assert all('error' not in result for result in results), results
    def identities(result):
        return {(r['sample'], r['group'], r['modality'], r['chemistry'], r['oligo_index'], r['sb_bc'])
                for sample in result['samples'] for r in sample['cell_identity_records']}
    assert identities(results[0]) == identities(results[1])
    assert identities(results[2]) == identities(results[3])
    assert results[0]['samples'][1]['split_targets']['VTD11_VTR12_K422_B_H3K27ac']['group'] == 'K422_B'
    records = results[0]['samples'][1]['cell_identity_records']
    assert any(r['oligo_index']=='07' and r['sb_bc']=='GAC' for r in records)
    assert results[0]['samples'][0]['sb_lookup_sha256'] == results[0]['samples'][1]['sb_lookup_sha256']
    assert results[-1]['samples'][0]['samplesheet_groups'] == ['RNA_only']
    assert results[-1]['samples'][1]['samplesheet_groups'] == ['DNA_only']


def test_authoritative_03_and_11_resolve_through_both_chemistries(tmp_path):
    results = parse_cases(tmp_path, [sheet({'B': group(sb_oligo_indices=['03', '11'])}, dna=mode)
                                     for mode in ('single', 'dual')])
    for result in results:
        assert 'error' not in result, result
        for row in result['samples']:
            expected = {'03':'TCA', '11':'ATG'} if row.get('dna_tagmentation') == 'dual' else {'03':'GACT', '11':'CTGA'}
            assert {r['oligo_index']:r['sb_bc'] for r in row['cell_identity_records']} == expected


def test_all_invalid_index_inputs_and_sequence_combinations(tmp_path):
    invalid = [True, False, 1.0, 0.1, 0, -1, '0', '-2', '+1', '1.0', ' 1', '1 ', '', 'one', '１２', None]
    cases = [sheet({'B':group(sb_oligo_indices=[entry])}) for entry in invalid]
    cases += [sheet({'B':group(sb_oligo_indices=entry)}) for entry in ([], '01', [1,'01'], ['13'])]
    for field in ('sb_barcodes', 'rna_sb_barcodes', 'dna_sb_barcodes'):
        cases.append(sheet({'B':group(sb_oligo_indices=['01'], **{field:['GCAT']})}))
    cases += [
        sheet({'B':group(rna_sb_barcodes=['GCAT'], dna_sb_barcodes=['AGT'])}),
        sheet({'B':group(rna_sb_barcodes=['AAAA'], dna_sb_barcodes=['GAT'])}),
        sheet({'B':group(rna_sb_barcodes=['GCAT'], dna_sb_barcodes=['AAA'])}),
        sheet({'B':group(sb_barcodes=['GCAT'])}), # legacy dual must be explicit
        sheet({'B':{'sb_oligo_indices':[1]}}), # DNA mark required
        sheet({'B':group(sb_oligo_indices=[1]), 'K422_B':group(sb_oligo_indices=['01'])}),
        sheet({'B':{'rna_sb_barcodes':['GCAT']}, 'K422_B':group(dna_sb_barcodes=['GAT'])}),
    ]
    results = parse_cases(tmp_path, cases)
    assert all('error' in result for result in results), results
    assert 'index sets differ' in results[-7]['error']


def test_namespaces_and_output_collisions_are_rejected(tmp_path):
    a = sheet({'C':group(sb_oligo_indices=[1])}, sample='A_B')
    b = sheet({'B_C':group(sb_oligo_indices=[2])}, sample='A')
    a['samples'].update(b['samples'])
    c = sheet({'B_C':{'dna_sb_barcodes':['GAT'], 'mark_barcodes':{'D':'GCCTCTAT'}},
               'B':{'dna_sb_barcodes':['AGT'], 'mark_barcodes':{'C_D':'GCCTCTAT'}}}, rna=False, sample='A')
    results = parse_cases(tmp_path, [a,c])
    assert 'namespace collision' in results[0]['error']
    assert 'output name collision' in results[1]['error']


def test_lookup_validation_and_extensibility(tmp_path):
    original = LOOKUP.read_text()
    valid = tmp_path / 'extended.tsv'
    valid.write_text(original + '13\tAAAA\tAAA\n100\tCCCC\tCCC\n')
    result = parse_cases(tmp_path/'valid', [sheet({'K422_B':group(sb_oligo_indices=[13,'00100'])})], valid)[0]
    assert 'error' not in result, result
    assert {r['oligo_index'] for r in result['samples'][0]['cell_identity_records']} == {'13','100'}
    mutations = [
        original + '001\tAAAA\tAAA\n', original + '13\tGCAT\tAAA\n', original + '13\tAAAA\tGAT\n',
        original + '13\t\tAAA\n', original + '13\taaaa\tAAA\n', original + '13\tAAAN\tAAA\n',
        original + '13\tAAA\tAAA\n', original + '13\tAAAA\tAAAA\n',
        original + '13\tAAAA\tAAA\textra\n', original.replace('oligo_index','index',1),
        original + '\n', original.splitlines()[0]+'\n',
    ]
    defaults = GROOVY_PARSE[GROOVY_PARSE.index('def defaults = '):GROOVY_PARSE.index('def contract = ')]
    files = []
    for i,text in enumerate(mutations):
        p=tmp_path/f'bad{i}.tsv';p.write_text(text);files.append(str(p))
    (tmp_path/'lookups.json').write_text(json.dumps(files))
    script='''def parser=new GroovyClassLoader().parseClass(new File(args[0]));
    def rows=new groovy.json.JsonSlurper().parse(new File(args[1])).collect { file ->
      try { parser.loadOligoLookup(new File(file)); return 'accepted' }
      catch(IllegalArgumentException e) { return e.message }
    }; println groovy.json.JsonOutput.toJson(rows)'''
    errors=json.loads(run(['java','-cp',find_nextflow_jar(),'groovy.ui.GroovyMain','-e',script,'--',
                          REPO/'lib/SamplesheetParser.groovy',tmp_path/'lookups.json']))
    assert 'accepted' not in errors


def test_exact_tonsil_tags_and_cross_chemistry_identity():
    dna = canonicalize_dna_fastq_comment('VTD11_VTR12','B','I:R:F:1:2:3:4:U',
       'CB:Z:TGCATATGCCAGTCATTGAGTAGTGACT\tSB:Z:TGCA\tMO:Z:GCCTCTAT','10')
    assert 'CB:Z:VTD11_VTR12_B_10_TATGCCAGTCATTGAGTAGTGACT' in dna
    assert 'XI:Z:VTD11_VTR12_B_10_TATGCCAGTCATTGAGTAGTGACT' in dna
    assert 'SB:Z:TGCA' in dna
    rna=canonicalize_fastq_comment('VTD11_VTR12','B',
       'CB:Z:GCTATTGACTCTTGATACGTTCCTCAAT\tSB:Z:GCTA\tUM:Z:ACGTACGTAC','07')
    assert 'CB:Z:VTD11_VTR12_B_07_TTGACTCTTGATACGTTCCTCAAT' in rna
    rna=canonicalize_fastq_comment('sample','K422_A',f'CB:Z:GCTA{LIGATIONS}\tSB:Z:GCTA','07')
    dna=canonicalize_dna_fastq_comment('sample','K422_A','I:R:F:1:2:3:4:U',f'CB:Z:GAC{LIGATIONS}\tSB:Z:GAC','07')
    assert next(t for t in rna.split() if t.startswith('CB:')) == next(t for t in dna.split() if t.startswith('CB:'))
    with pytest.raises(ValueError):
        canonicalize_fastq_comment('sample','A',f'CB:Z:GCTA{LIGATIONS}\tSB:Z:AGCTA','07')


def split_fixture(root, modality, mode, sb_sequences=('GCAT','CGAT'), indices=('01','02'),
                  sample='VTD11_VTR12', sb_map=None, mo_map=None):
    root.mkdir(parents=True,exist_ok=True)
    records=[]
    for index, sb in zip(indices,sb_sequences):
        for read in range(4):
            qname=f'I:R:F:1:1:{int(index)*100+read}:1:U'
            comment=f'CB:Z:{sb}{LIGATIONS}\tSB:Z:{sb}\tUM:Z:ACGTACGTAC\tMO:Z:GCCTCTAT'
            records.append((qname,comment))
    records.append(('I:R:F:1:1:9999:1:U','CB:Z:NoMatch\tSB:Z:NoMatch\tUM:Z:ACGTACGTAC\tMO:Z:GCCTCTAT'))
    for role in ('r1','r2'):
        (root/f'{role}.fastq').write_text(''.join(f'@{q} {c}\nACGTACGT\n+\nIIIIIIII\n' for q,c in records))
    (root/'sb.tsv').write_text('sample\tsb_group\tsb_bc\toligo_index\n'+''.join(
        f'{sample}\tK422_B\t{sb}\t{index}\n' for sb,index in zip(sb_sequences,indices)))
    (root/'mo.tsv').write_text(f'sample\tsb_group\tmark\tmo_bc\n{sample}\tK422_B\tH3K27ac\tGCCTCTAT\n')
    if sb_map:
        shutil.copyfile(sb_map, root/'sb.tsv')
    if mo_map:
        shutil.copyfile(mo_map, root/'mo.tsv')
    command=[sys.executable,REPO/f'bin/run_split_reads_{modality}.py','--mode',mode,'--script',
             REPO/'scripts/core_runtime/Split_ReadsV2.codon','--r1',root/'r1.fastq','--r2',root/'r2.fastq',
             '--sb-group-map',root/'sb.tsv','--sample',sample,'--library-name','TEST','--output-dir',root/'out']
    if modality=='dna':command+=['--mo-map',root/'mo.tsv']
    run(command, env=dict(os.environ, TMPDIR=str(root)))
    return root/'out'


def test_sample_named_sample_uses_real_derived_headers(tmp_path):
    parsed = parse_cases(tmp_path/'contract', [sheet({'K422_B':group(sb_oligo_indices=[1,2])}, sample='sample')])[0]
    assert 'error' not in parsed, parsed
    for row in parsed['samples']:
        modality = row['modality']
        sequences = ('GCAT','CGAT') if modality == 'rna' else ('GAT','AGT')
        outputs = {mode:split_fixture(tmp_path/modality/mode, modality, mode, sequences, sample='sample',
                                      sb_map=row['sb_group_map'], mo_map=row.get('mo_map'))
                   for mode in ('mock','real')}
        for path in outputs['real'].glob('*.fastq'):
            assert path.read_text() == (outputs['mock']/path.name).read_text()
            assert 'CB:Z:sample_K422_B_01_' in path.read_text()
            assert 'CB:Z:sample_K422_B_02_' in path.read_text()
        assert not any('sb_group' in path.name for path in outputs['mock'].iterdir())


@pytest.mark.parametrize('modality,tag', [('rna','SB'), ('dna','SB'), ('dna','MO'), ('rna','UM')])
def test_real_mock_split_rejects_mate_tag_mismatches(tmp_path, modality, tag):
    root = tmp_path/modality
    sequences = ('GCAT','CGAT') if modality == 'rna' else ('GAT','AGT')
    split_fixture(root, modality, 'mock', sequences)
    lines = (root/'r2.fastq').read_text().splitlines()
    if tag == 'SB':
        # Both raw SB and raw CB change, keeping ligations identical. Checking
        # only rewritten tags would hide this mismatch by copying R1's index.
        lines[0] = lines[0].replace(f'CB:Z:{sequences[0]}', f'CB:Z:{sequences[1]}').replace(
            f'SB:Z:{sequences[0]}', f'SB:Z:{sequences[1]}')
    else:
        lines[0] = lines[0].replace('MO:Z:GCCTCTAT','MO:Z:AGGCTATA') if tag == 'MO' else lines[0].replace(
            'UM:Z:ACGTACGTAC','UM:Z:TGCATGCATG')
    (root/'r2.fastq').write_text('\n'.join(lines)+'\n')
    for mode in ('mock','real'):
        command = [sys.executable, REPO/f'bin/run_split_reads_{modality}.py', '--mode',mode,
                   '--script',REPO/'scripts/core_runtime/Split_ReadsV2.codon', '--sample','VTD11_VTR12',
                   '--library-name','TEST', '--r1',root/'r1.fastq', '--r2',root/'r2.fastq',
                   '--sb-group-map',root/'sb.tsv', '--output-dir',root/mode]
        if modality == 'dna':
            command += ['--mo-map',root/'mo.tsv']
        result = subprocess.run(list(map(str,command)), capture_output=True, text=True,
                                env=dict(os.environ, TMPDIR=str(root)))
        assert result.returncode != 0, (modality,tag,mode)
        assert 'Mate' in result.stdout+result.stderr


def test_real_mock_split_fq_sam_parity_and_nomatch(tmp_path):
    assert shutil.which('codon'), 'Real Codon required'
    for modality,sequences in [('rna',('GCAT','CGAT')),('dna',('GAT','AGT'))]:
        output={mode:split_fixture(tmp_path/modality/mode,modality,mode,sequences,('13','100')) for mode in ('mock','real')}
        for path in output['real'].glob('*.fastq'):
            assert path.read_text()==(output['mock']/path.name).read_text()
            headers=path.read_text().splitlines()[::4]
            assert len(headers)==8 and all('NoMatch' not in h for h in headers)
            assert {h.split('CB:Z:',1)[1].split()[0] for h in headers}=={
                f'VTD11_VTR12_K422_B_{idx}_{LIGATIONS}' for idx in ('13','100')}
        if modality=='rna':
            for mode in ('mock','real'):
                out=output[mode]
                run([sys.executable,REPO/'bin/run_fq_to_sam.py','--mode',mode,'--script',REPO/'scripts/core_runtime/FqToSAM.codon',
                     '--r1',out/'VTD11_VTR12_K422_B_R1.fastq','--r2',out/'VTD11_VTR12_K422_B_R2.fastq','--output-sam',out/'rna.sam'])
            assert (output['mock']/'rna.sam').read_text()==(output['real']/'rna.sam').read_text()
            text=(output['real']/'rna.sam').read_text()
            assert 'CR:Z:' not in text and 'UR:Z:ACGTACGTAC' in text


def test_mock_solo_summary_matches_the_full_identity_matrix(tmp_path):
    out = split_fixture(tmp_path/'split', 'rna', 'mock')
    usam = tmp_path/'rna.sam'
    run([sys.executable,REPO/'bin/run_fq_to_sam.py','--mode','mock','--script',REPO/'scripts/core_runtime/FqToSAM.codon',
         '--r1',out/'VTD11_VTR12_K422_B_R1.fastq','--r2',out/'VTD11_VTR12_K422_B_R2.fastq','--output-sam',usam])
    solo = tmp_path/'Solo.outGeneFull'
    run([sys.executable, REPO/'bin/mock_starsolo_outputs.py', usam, solo])
    summary = dict(csv.reader((solo/'Summary.csv').open()))
    barcodes = (solo/'filtered/barcodes.tsv').read_text().splitlines()
    assert summary['Estimated Number of Cells'] == str(len(barcodes)) == '2'
    assert summary['UMIs in Cells'] == '2' and summary['Number of Reads'] == '8'
    assert float(summary['Sequencing Saturation']) == 0.75
    matrix = [line for line in (solo/'filtered/matrix.mtx').read_text().splitlines() if not line.startswith('%')]
    assert matrix == ['1 2 2', '1 1 1', '1 2 1']
    usam.write_text('@HD\tVN:1.6\tSO:unsorted\n')
    run([sys.executable, REPO/'bin/mock_starsolo_outputs.py', usam, solo])
    summary = dict(csv.reader((solo/'Summary.csv').open()))
    assert summary['Estimated Number of Cells'] == summary['UMIs in Cells'] == '0'
    assert (solo/'filtered/barcodes.tsv').read_bytes() == b''


@pytest.mark.parametrize('sample,group_name', [
    ('NoMatchControl', 'K422_B'), ('NoMatch', 'K422_B'),
    ('sample_id', 'NoMatchControl'), ('sample_id', 'NoMatch'),
])
def test_real_mock_fq_to_sam_keeps_nomatch_names_and_excludes_failed_tags(tmp_path, sample, group_name):
    assert shutil.which('codon'), 'Real Codon required'
    full = f'{sample}_{group_name}_01_{LIGATIONS}'
    comment = f'CB:Z:{full}\tXI:Z:{full}\tSB:Z:GCAT\tUM:Z:ACGTACGTAC\tMO:Z:GCCTCTAT\tRG:Z:{full}'
    mates = {1:[('valid', comment)], 2:[('valid', comment)]}
    # A failed barcode on either mate must still exclude the entire pair.
    for tag, value in [('CB', full), ('SB', 'GCAT'), ('MO', 'GCCTCTAT')]:
        for failed_mate in (1, 2):
            qname = f'failed_{tag}_{failed_mate}'
            for mate in (1, 2):
                tags = comment.replace(f'{tag}:Z:{value}', f'{tag}:Z:NoMatch') if mate == failed_mate else comment
                mates[mate].append((qname, tags))
    for mate, records in mates.items():
        (tmp_path/f'R{mate}.fastq').write_text(''.join(
            f'@{name} {tags}\nACGTACGT\n+\nIIIIIIII\n' for name, tags in records))
    for mode in ('mock', 'real'):
        run([sys.executable, REPO/'bin/run_fq_to_sam.py', '--mode',mode,
             '--script',REPO/'scripts/core_runtime/FqToSAM.codon',
             '--r1',tmp_path/'R1.fastq', '--r2',tmp_path/'R2.fastq', '--output-sam',tmp_path/f'{mode}.sam'])
        records = [line.split('\t') for line in (tmp_path/f'{mode}.sam').read_text().splitlines()
                   if not line.startswith('@')]
        assert len(records) == 2, (sample, group_name, mode)
        assert {record[0] for record in records} == {'valid'}
        for record in records:
            tags = dict(token.split(':',2)[::2] for token in record[11:])
            assert tags['CB'] == tags['XI'] == full
            assert tags['UR'] == tags['UM'] == 'ACGTACGTAC'
    assert (tmp_path/'mock.sam').read_bytes() == (tmp_path/'real.sam').read_bytes()


def test_real_star_rejected_umis_and_unassigned_reads_keep_identity(tmp_path):
    assert shutil.which('STAR') and shutil.which('samtools') and shutil.which('codon'), 'Real RNA tools required'
    genome = ''.join(random.Random(89).choices('ACGT', k=100000))
    (tmp_path/'genome.fa').write_text('>chr1\n'+genome+'\n')
    # g2 and g3 overlap exactly, giving a unique alignment with an ambiguous
    # gene assignment. The outside case maps beyond all annotated genes.
    (tmp_path/'genes.gtf').write_text(''.join(
        f'chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene}"; transcript_id "t{gene}"; gene_name "{gene}";\n'
        for gene,start,end in [('g1',101,2000),('g2',4001,6000),('g3',4001,6000)]))
    index = tmp_path/'index'; index.mkdir()
    run(['STAR','--runMode','genomeGenerate','--genomeDir',index,'--genomeFastaFiles',tmp_path/'genome.fa',
         '--sjdbGTFfile',tmp_path/'genes.gtf','--sjdbOverhang','49','--genomeSAindexNbases','5',
         '--genomeChrBinNbits','10','--runThreadN','1'], cwd=tmp_path)
    cases = [('valid','ACGTACGTAC',200), ('N','ACGTNCGTAC',200)]
    cases += [(f'poly{base}',base*10,200) for base in 'ACGT']
    cases += [('outside','TGCATGCATG',8000), ('ambiguous','GATCGATCGA',4500)]
    reads = {}
    for idx, sb in [('01','GCAT'), ('02','CGAT')]:
        for kind, umi, start in cases:
            reads[f'idx{idx}_{kind}'] = {'cell':f'RNA_sample_K422_A_{idx}_{LIGATIONS}',
                                       'umi':umi, 'sb':sb, 'kind':kind, 'start':start}
    complement = str.maketrans('ACGT','TGCA')
    for mate, offset in [(1,0), (2,100)]:
        records = []
        for qname, read in reads.items():
            sequence = genome[read['start']+offset:read['start']+offset+50]
            if mate == 2:
                sequence = sequence.translate(complement)[::-1]
            comment = f"CB:Z:{read['cell']}\tXI:Z:{read['cell']}\tSB:Z:{read['sb']}\tUM:Z:{read['umi']}\tRG:Z:{read['cell']}"
            records.append(f'@{qname} {comment}\n{sequence}\n+\n'+ 'I'*50+'\n')
        (tmp_path/f'R{mate}.fastq').write_text(''.join(records))
    usam = tmp_path/'rna.sam'
    run([sys.executable,REPO/'bin/run_fq_to_sam.py','--mode','real',
         '--script',REPO/'scripts/core_runtime/FqToSAM.codon','--r1',tmp_path/'R1.fastq',
         '--r2',tmp_path/'R2.fastq','--output-sam',usam])
    script = (REPO/'scripts/core_runtime/RNA_STARSOLO_ALIGN.sh').read_text().replace(
        '--soloCellFilter EmptyDrops_CR 15000 0.99 30 20000 90000 50 0.01 20000 0.05 10000',
        '--soloCellFilter TopCells 2')
    # Capture STAR's actual tags before the production normalizer to verify
    # that every retained alignment and UB value survives unchanged.
    script = script.replace('  | "${PYTHON3_BIN}"',
                            '  | tee "${outdir}/star_before_normalization.sam" \\\n  | "${PYTHON3_BIN}"')
    (tmp_path/'align.sh').write_text(script)
    shutil.copyfile(REPO/'scripts/core_runtime/NormalizeRnaBamTags.py',tmp_path/'NormalizeRnaBamTags.py')
    run(['bash',tmp_path/'align.sh','fixture',usam,index,tmp_path,'1'], cwd=tmp_path)
    raw = [line.split('\t') for line in (tmp_path/'star_before_normalization.sam').read_text().splitlines()
           if not line.startswith('@')]
    normalized = [line.split('\t') for line in run(
        ['samtools','view',tmp_path/'fixture.Aligned.sortedByCoord.out.bam']).splitlines()]
    assert len(raw) == len(normalized) == len(reads)*2
    assert {r[0] for r in normalized} == set(reads)
    evidence = []
    for before, after in zip(raw, normalized):
        assert before[:11] == after[:11]
        assert [t for t in before[11:] if not t.startswith(('CB:', 'XI:'))] == [
            t for t in after[11:] if not t.startswith(('CB:', 'XI:'))]
        read = reads[after[0]]
        raw_cb = [t[5:] for t in before[11:] if t.startswith('CB:Z:')]
        tags = dict(t.split(':',2)[::2] for t in after[11:])
        assert tags['CB'] == tags['XI'] == read['cell']
        assert sum(t.startswith('CB:') for t in after[11:]) == sum(t.startswith('XI:') for t in after[11:]) == 1
        assert tags['UR'] == tags['UM'] == read['umi']
        if read['kind'] == 'N' or read['kind'].startswith('poly'):
            assert raw_cb == [read['cell'], '-']
            assert tags['UB'] == '-'
        elif read['kind'] == 'valid':
            assert tags['UB'] == read['umi']
        elif read['kind'] in ('outside','ambiguous'):
            assert tags['GX'] == '-'
            assert tags['UB'] == read['umi']
        evidence.append({'qname':after[0], 'flag':int(after[1]), 'kind':read['kind'],
                         'raw_CB':raw_cb, 'CB':tags['CB'], 'XI':tags['XI'],
                         'UM':tags['UM'], 'UB':tags['UB'], 'GX':tags.get('GX')})
    solo = tmp_path/'fixture.Solo.outGeneFull'
    expected = [f'RNA_sample_K422_A_{idx}_{LIGATIONS}' for idx in ('01','02')]
    matrices = {}
    for layer in ('raw','filtered'):
        assert (solo/layer/'barcodes.tsv').read_text().splitlines() == expected
        features = [row.split('\t')[0] for row in (solo/layer/'features.tsv').read_text().splitlines()]
        matrix = [line for line in (solo/layer/'matrix.mtx').read_text().splitlines() if not line.startswith('%')]
        assert matrix[0] == '3 2 2'
        assert [(features[int(g)-1],int(c),int(n)) for g,c,n in (line.split() for line in matrix[1:])] == [
            ('g1',1,1), ('g1',2,1)]
        matrices[layer] = matrix
    em_matrix = [line for line in (solo/'raw/UniqueAndMult-EM.mtx').read_text().splitlines()
                 if not line.startswith('%')]
    assert em_matrix[0] == '3 2 6'
    assert {(features[int(g)-1],int(c)):float(n) for g,c,n in (line.split() for line in em_matrix[1:])} == {
        (gene,cell):count for cell in (1,2) for gene,count in [('g1',1.0),('g2',0.5),('g3',0.5)]}
    matrices['raw_EM'] = em_matrix
    stats = list(csv.DictReader((solo/'CellReads.stats').open(),delimiter='\t'))
    assert {r['CB']:int(r['nUMIunique']) for r in stats if r['CB'] != 'CBnotInPasslist'} == dict.fromkeys(expected,1)
    # A second non-placeholder cell identity must still fail the executable
    # normalizer, even on a STAR record that also contains its '-' placeholder.
    rejected = next(r for r in raw if r[0] == 'idx01_N')
    corrupt = '\t'.join(rejected + [f'CB:Z:RNA_sample_K422_A_02_{LIGATIONS}'])+'\n'
    conflict = subprocess.run([sys.executable,REPO/'scripts/core_runtime/NormalizeRnaBamTags.py'],
                              input=corrupt,capture_output=True,text=True)
    assert conflict.returncode != 0 and 'Conflicting CB' in conflict.stderr
    # Called-cell BAM selection is independent of whether STAR accepted the
    # UMI. Keep all existing proper-pair alignments for these two called cells.
    (tmp_path/'canonical.txt').write_text('chr1\n')
    run(['bash',REPO/'scripts/core_runtime/RNA_FILTERED_BAM.sh','fixture',solo,
         tmp_path/'fixture.Aligned.sortedByCoord.out.bam',tmp_path/'canonical.txt',tmp_path,'1'])
    filtered = [line.split('\t') for line in run(
        ['samtools','view',tmp_path/'fixture.filtered_cells.bam']).splitlines()]
    assert filtered == normalized
    retention = list(csv.DictReader((tmp_path/'fixture.rna_filter_retention.tsv').open(),delimiter='\t'))
    assert next(int(r['pairs']) for r in retention if r['metric'] == 'called_cell_pairs') == len(reads)
    (tmp_path/'evidence.json').write_text(json.dumps(
        {'alignments':evidence,'matrices':matrices,'cell_stats':stats,'retention':retention,
         'filtered_alignment_count':len(filtered),'genuine_conflict_rejected':True},indent=2)+'\n')


def test_real_star_two_indices_independent_counts_and_filtered_bam(tmp_path):
    assert shutil.which('STAR') and shutil.which('samtools'), 'Real STAR/samtools required'
    genome=''.join(random.Random(81).choices('ACGT',k=100000))
    (tmp_path/'genome.fa').write_text('>chr1\n'+genome+'\n')
    (tmp_path/'genes.gtf').write_text('chr1\ttest\texon\t101\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "G1";\n')
    index=tmp_path/'index';index.mkdir()
    run(['STAR','--runMode','genomeGenerate','--genomeDir',index,'--genomeFastaFiles',tmp_path/'genome.fa',
         '--sjdbGTFfile',tmp_path/'genes.gtf','--sjdbOverhang','49','--genomeSAindexNbases','5','--genomeChrBinNbits','10','--runThreadN','1'],cwd=tmp_path)
    # Use real splitter + real FqToSAM, with two different SBs and the same UMI.
    out=split_fixture(tmp_path/'split','rna','real')
    complement=str.maketrans('ACGT','TGCA')
    for mate,start,reverse in [('R1',200,False),('R2',300,True)]:
        p=out/f'VTD11_VTR12_K422_B_{mate}.fastq';lines=p.read_text().splitlines()
        sequence=genome[start:start+50]
        if reverse:sequence=sequence.translate(complement)[::-1]
        for i in range(0,len(lines),4):lines[i+1]=sequence;lines[i+3]='I'*50
        p.write_text('\n'.join(lines)+'\n')
    usam=tmp_path/'rna.sam'
    run([sys.executable,REPO/'bin/run_fq_to_sam.py','--mode','real','--script',REPO/'scripts/core_runtime/FqToSAM.codon',
         '--r1',out/'VTD11_VTR12_K422_B_R1.fastq','--r2',out/'VTD11_VTR12_K422_B_R2.fastq','--output-sam',usam])
    script=(REPO/'scripts/core_runtime/RNA_STARSOLO_ALIGN.sh').read_text().replace(
       '--soloCellFilter EmptyDrops_CR 15000 0.99 30 20000 90000 50 0.01 20000 0.05 10000','--soloCellFilter TopCells 2')
    (tmp_path/'align.sh').write_text(script)
    shutil.copyfile(REPO/'scripts/core_runtime/NormalizeRnaBamTags.py',tmp_path/'NormalizeRnaBamTags.py')
    run(['bash',tmp_path/'align.sh','fixture',usam,index,tmp_path,'1'],cwd=tmp_path)
    expected=[f'VTD11_VTR12_K422_B_{idx}_{LIGATIONS}' for idx in ('01','02')]
    solo=tmp_path/'fixture.Solo.outGeneFull'
    for layer in ('raw','filtered'):
        assert (solo/layer/'barcodes.tsv').read_text().splitlines()==expected
        matrix=[l for l in (solo/layer/'matrix.mtx').read_text().splitlines() if not l.startswith('%')]
        assert matrix==['1 2 2','1 1 1','1 2 1']
    stats=list(csv.DictReader((solo/'CellReads.stats').open(),delimiter='\t'))
    assert {r['CB'] for r in stats if r['CB']!='CBnotInPasslist'}==set(expected)
    assert [int(r['nUMIunique']) for r in stats if r['CB']!='CBnotInPasslist']==[1,1]
    bam=tmp_path/'fixture.Aligned.sortedByCoord.out.bam'
    records=run(['samtools','view',bam]).splitlines()
    assert len(records)==16
    for record in records:
        tokens=record.split('\t')[11:]
        tags=dict(t.split(':',2)[::2] for t in tokens)
        assert tags['CB']==tags['XI'] and tags['CB'] in expected
        assert sum(t.startswith('CB:') for t in tokens)==1
        assert tags['UR']==tags['UB']=='ACGTACGTAC'
    # Production called-cell filter and its audit use exactly the same XI set.
    (tmp_path/'canonical.txt').write_text('chr1\n')
    run(['bash',REPO/'scripts/core_runtime/RNA_FILTERED_BAM.sh','fixture',solo,bam,tmp_path/'canonical.txt',tmp_path,'1'])
    filtered=run(['samtools','view',tmp_path/'fixture.filtered_cells.bam']).splitlines()
    assert len(filtered)==16
    retention=list(csv.DictReader((tmp_path/'fixture.rna_filter_retention.tsv').open(),delimiter='\t'))
    assert next(int(r['pairs']) for r in retention if r['metric']=='called_cell_pairs')==8


def test_identity_normalization_preserves_secondary_and_supplementary():
    spec=importlib.util.spec_from_file_location('normalize',REPO/'scripts/core_runtime/NormalizeRnaBamTags.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    full=f'sample_K422_A_100_{LIGATIONS}'
    for flag in (99,147,355,403,2147,2195):
        for cells in (f'CB:Z:{full}\tXI:Z:{full}\tCB:Z:{full}',
                      f'CB:Z:{full}\tXI:Z:{full}\tCB:Z:-',
                      f'CB:Z:-\tCB:Z:{full}\tXI:Z:{full}'):
            line=f'q\t{flag}\tchr1\t1\t60\t4M\t=\t5\t0\tACGT\tIIII\t{cells}\tSB:Z:GCAT\tUR:Z:ACGTACGTAC\tUB:Z:-\n'
            normalized=module.normalize_record(line)
            assert f'\t{flag}\t' in normalized
            assert normalized.count('CB:Z:')==normalized.count('XI:Z:')==1
            assert f'CB:Z:{full}' in normalized and f'XI:Z:{full}' in normalized
            assert 'UB:Z:-' in normalized


def test_identity_normalization_rejects_nonplaceholder_conflicts_and_missing_input_cb():
    spec=importlib.util.spec_from_file_location('normalize',REPO/'scripts/core_runtime/NormalizeRnaBamTags.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    full = f'sample_K422_A_01_{LIGATIONS}'
    other = f'sample_K422_A_02_{LIGATIONS}'
    for cells in (f'CB:Z:{full}\tXI:Z:{full}\tCB:Z:-\tCB:Z:{other}',
                  f'CB:Z:{full}\tXI:Z:{full}\tXI:Z:{other}\tCB:Z:-',
                  f'CB:Z:{full}\tXI:Z:{other}\tCB:Z:-',
                  f'CB:Z:-\tXI:Z:{full}', f'CB:Z:{full}\tXI:Z:-'):
        line=f'q\t99\tchr1\t1\t60\t4M\t=\t5\t0\tACGT\tIIII\t{cells}\n'
        with pytest.raises(ValueError):
            module.normalize_record(line)


def test_real_dna_alignment_preserves_full_tags_on_both_mates(tmp_path):
    assert shutil.which('bwa-mem2') and shutil.which('samtools'), 'Real bwa-mem2/samtools required'
    genome=''.join(random.Random(91).choices('ACGT',k=100000))
    fasta=tmp_path/'genome.fa';fasta.write_text('>chr1\n'+genome+'\n')
    run(['bwa-mem2','index',fasta],cwd=tmp_path)
    out=split_fixture(tmp_path/'split','dna','real',('GAT','AGT'))
    complement=str.maketrans('ACGT','TGCA')
    for mate,start,reverse in [('R1',200,False),('R2',400,True)]:
        p=out/f'VTD11_VTR12_K422_B_H3K27ac_{mate}.fastq';lines=p.read_text().splitlines()
        sequence=genome[start:start+100]
        if reverse:sequence=sequence.translate(complement)[::-1]
        for i in range(0,len(lines),4):lines[i+1]=sequence;lines[i+3]='I'*100
        # Sixteen pairs provide enough insert observations for bwa-mem2's
        # ordinary proper-pair inference, keeping production options intact.
        extra=lines.copy()
        for i in range(0,len(extra),4):
            extra[i]=extra[i].replace(':1:U ',':2:U ',1)
        p.write_text('\n'.join(lines+extra)+'\n')
    blacklist=tmp_path/'blacklist.bed';blacklist.write_text('chr1\t90000\t90001\n')
    environment=dict(os.environ,ALIGN_DNA_THREADS='1',ALIGN_DNA_VIEW_THREADS='1',ALIGN_DNA_SORT_THREADS='1',ALIGN_DNA_SORT_MEM='64M')
    run(['bash',REPO/'scripts/core_runtime/AlignDNA.sh','H3K27ac','VTD11_VTR12_K422_B',
         out/'VTD11_VTR12_K422_B_H3K27ac_R1.fastq',out/'VTD11_VTR12_K422_B_H3K27ac_R2.fastq',blacklist,
         out/'SAM_RG_Header_VTD11_VTR12_K422_B_H3K27ac.tsv',fasta,'100000',tmp_path],cwd=tmp_path,env=environment)
    records=run(['samtools','view',tmp_path/'VTD11_VTR12_K422_B_H3K27ac.bam']).splitlines()
    assert len(records)==32
    expected={f'VTD11_VTR12_K422_B_{index}_{LIGATIONS}' for index in ('01','02')}
    for record in records:
        tags=dict(t.split(':',2)[::2] for t in record.split('\t')[11:])
        assert tags['CB']==tags['XI'] and tags['CB'] in expected
        assert tags['SB'] in ('GAT','AGT') and tags['RG']=='I:R:F:L1'
    header=run(['samtools','view','-H',tmp_path/'VTD11_VTR12_K422_B_H3K27ac.bam'])
    assert '\tLB:TEST\tPU:I:R:F:L1' in header
