#!/usr/bin/env python3
"""Create a bounded synchronized development samplesheet and raw FASTQ subset.

Usage: make_cell_identity_subset.py --samplesheet ORIGINAL --output-dir NEW
       [--sample ID] [--pairs-per-read-set 1000]
Requires PyYAML in the development environment, never modifies the source.
"""
import argparse
import gzip
from itertools import zip_longest
from pathlib import Path

import yaml

from tresflow_fastq_utils import fastq_iter, normalize_qname, parse_header


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--sample', action='append', help='Restrict to these sample keys')
    parser.add_argument('--pairs-per-read-set', type=int, default=1000)
    args = parser.parse_args()
    if args.pairs_per_read_set <= 0:
        parser.error('--pairs-per-read-set must be positive')
    # An exclusive root prevents overwriting a prior development run.
    args.output_dir.mkdir(parents=True, exist_ok=False)
    source = args.samplesheet.resolve()
    sheet = yaml.safe_load(source.read_text())
    sheet['runtime']['tmpdir'] = str((args.output_dir / 'runtime_tmp').resolve())
    if args.sample:
        missing = set(args.sample) - sheet['samples'].keys()
        if missing:
            parser.error(f'Unknown samples: {sorted(missing)}')
        sheet['samples'] = {s: sheet['samples'][s] for s in args.sample}
    for key, value in sheet['references'].items():
        if key not in ('species', 'dna_effective_genome_size') and isinstance(value, str):
            sheet['references'][key] = str((source.parent / value).resolve())
    for sample, config in sheet['samples'].items():
        for modality in ('rna', 'dna'):
            if modality not in config:
                continue
            reads = config[modality]['reads']
            paths = {role: value if isinstance(value, list) else value.split(',') for role, value in reads.items()}
            sizes = {len(values) for values in paths.values()}
            if len(sizes) != 1:
                raise ValueError(f'Conflicting read-set counts for {sample}/{modality}')
            outputs = {role: [] for role in paths}
            for i in range(next(iter(sizes))):
                roles = list(paths)
                inputs = [fastq_iter((source.parent / paths[role][i].strip()).resolve()) for role in roles]
                directory = args.output_dir / sample / modality / str(i+1)
                directory.mkdir(parents=True)
                handles = []
                try:
                    for role in roles:
                        path = directory / f'{role}.fastq.gz'
                        handles.append(gzip.open(path, 'wt'))
                        outputs[role].append(str(path.resolve()))
                    for count, records in enumerate(zip_longest(*inputs)):
                        if any(r is None for r in records):
                            raise ValueError(f'Unequal EOF in {sample}/{modality}/{i+1}')
                        names = {normalize_qname(parse_header(r[0])[0]) for r in records}
                        if len(names) != 1:
                            raise ValueError(f'QNAME mismatch in {sample}/{modality}/{i+1}/{count+1}')
                        for handle, record in zip(handles, records):
                            handle.write('\n'.join(record)+'\n')
                        if count+1 == args.pairs_per_read_set:
                            break
                finally:
                    for handle in handles:
                        handle.close()
                    for input_stream in inputs:
                        input_stream.close()
            config[modality]['reads'] = outputs
    (args.output_dir / 'samplesheet.dev.yaml').write_text(yaml.safe_dump(sheet, sort_keys=False))


if __name__ == '__main__':
    main()
