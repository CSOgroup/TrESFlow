#!/usr/bin/env python3
"""Export routed RNA pairs for bulk consumers; invoked only in preprocessing mode."""
import argparse
import csv
import hashlib
import gzip
import json
import sys
import subprocess
from itertools import zip_longest
from functools import lru_cache
from pathlib import Path

from tresflow_fastq_utils import fastq_iter, normalize_qname, parse_header


REQUIRED_TAGS = ('CB', 'XI', 'SB', 'UM', 'RG')
# Keep both compressors fed: larger per-mate writes can block on one pipe
# long enough to starve the other, particularly with single-threaded pigz.
WRITE_BUFFER_SIZE = 1 << 14


def compressor_threads(cpus):
    """Reserve conversion plus both pigz workers and their read/write threads."""
    if cpus < 1:
        raise ValueError('RNA export CPU budget must be positive')
    if cpus < 3:
        return 0  # Serial in-process gzip avoids oversubscribing small tasks.
    # pigz -p 1 is single-threaded; -p N>1 also has two I/O threads.
    return max(1, (cpus - 1) // 2 - 2)


def parse_tags(comment):
    # Match find_tag_value's first-occurrence and malformed-token semantics,
    # without splitting/scanning the whole comment once for every tag.
    tags = {}
    for token in comment.split():
        fields = token.split(':', 2)
        key = fields[0]
        if key in REQUIRED_TAGS and key not in tags:
            tags[key] = fields[2] if len(fields) == 3 else ''
    return tags


@lru_cache(maxsize=65536)
def pool_identifier(cb, sb):
    # Fixed length, delimiter-safe, stable across technical chunks. The mapping
    # sidecar and original comments retain the full identity and physical SB.
    return 'P' + hashlib.sha256(json.dumps([cb, sb], separators=(',', ':')).encode()).hexdigest()


def convert_pair(r1, r2, prefix_length, umi_length, minimum_length):
    name1, comment1 = parse_header(r1[0])
    name2, comment2 = parse_header(r2[0])
    original = normalize_qname(name1)
    if original != normalize_qname(name2):
        raise ValueError(f'RNA export mate names differ: {name1} != {name2}')
    tags = parse_tags(comment1)
    mate_tags = tags if comment1 == comment2 else parse_tags(comment2)
    if any(not tags.get(tag) or tags[tag] == 'NoMatch' or tags[tag] != mate_tags.get(tag, '') for tag in REQUIRED_TAGS):
        raise ValueError(f'RNA export requires matching accepted CB/XI/SB/UM/RG tags: {original}')
    if tags['CB'] != tags['XI'] or len(tags['UM']) != umi_length or tags['UM'].strip('ACGTNacgtn'):
        raise ValueError(f'Invalid RNA export CB/XI or UM identity: {original}')
    pool = pool_identifier(tags['CB'], tags['SB'])
    # UM is already reverse-complemented by TAG_RNA_UMI. Copy it verbatim.
    name = f"{original}_{pool}_{tags['UM']}"
    sequence2, quality2 = r2[1][prefix_length:], r2[3][prefix_length:]
    if len(r1[1]) < minimum_length or len(sequence2) < minimum_length:
        return None, pool, tags
    return (
        f'@{name} {comment1}\n{r1[1]}\n+\n{r1[3]}\n',
        f'@{name} {comment2}\n{sequence2}\n+\n{quality2}\n',
    ), pool, tags


def export_fastqs(manifest, output_dir, sample, prefix_length, umi_length, minimum_length, pigz, threads):
    mate_threads = compressor_threads(threads)
    if mate_threads:
        print(f'RNA export CPU budget={threads}: converter=1, pigz -p {mate_threads} per mate (including I/O thread headroom)', file=sys.stderr)
    else:
        print(f'RNA export CPU budget={threads}: serial conversion/gzip in one process', file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics, pools = [], {}
    for split in manifest:
        paths = [output_dir / (Path(split[role]).name + '.gz') for role in ('r1', 'r2')]
        input_pairs = retained = 0
        processes, handles, streams = [], [], []
        buffers = [[], []]
        buffered_size = 0
        split_pools = {}

        def flush_buffers():
            for stream, records in zip(streams, buffers):
                if records:
                    stream.write(''.join(records))
                    records.clear()

        try:
            for path in paths:
                if mate_threads:
                    handles.append(path.open('wb'))
                    processes.append(subprocess.Popen([str(pigz), '-c', '-p', str(mate_threads)], stdin=subprocess.PIPE, stdout=handles[-1], text=True, bufsize=WRITE_BUFFER_SIZE))
                    streams.append(processes[-1].stdin)
                else:
                    streams.append(gzip.open(path, 'wt', compresslevel=6))
            for r1, r2 in zip_longest(fastq_iter(Path(split['r1'])), fastq_iter(Path(split['r2']))):
                if r1 is None or r2 is None:
                    raise ValueError('RNA export mates have unequal EOF')
                input_pairs += 1
                pair, pool, tags = convert_pair(r1, r2, prefix_length, umi_length, minimum_length)
                if pool not in split_pools:
                    split_pools[pool] = [sample, split['group'], split['split_name'], pool, tags['CB'], tags['XI'], tags['SB']]
                if pair is not None:
                    buffers[0].append(pair[0])
                    buffers[1].append(pair[1])
                    buffered_size += len(pair[0]) + len(pair[1])
                    if buffered_size >= 2 * WRITE_BUFFER_SIZE:
                        flush_buffers()
                        buffered_size = 0
                    retained += 1
            flush_buffers()
            for stream in streams:
                stream.close()
            codes = [process.wait() for process in processes]
            if any(codes):
                raise RuntimeError(f'RNA export pigz failed: {codes}')
        finally:
            failed = sys.exc_info()[0] is not None
            for process in processes:
                if process.poll() is None:
                    process.kill()
                    process.wait()
            for stream in streams:
                if not stream.closed:
                    try:
                        stream.close()
                    except OSError:
                        # A dead compressor can also fail while closing its
                        # pipe; retain the original conversion/I/O exception.
                        if not failed:
                            raise
            for handle in handles:
                handle.close()
        pools.update(((split['split_name'], pool), row) for pool, row in split_pools.items())
        for metric, count in [('export_input_pairs', input_pairs), ('exported_pairs', retained), ('export_short_pairs', input_pairs - retained)]:
            metrics.append([sample, 'rna', split['group'], 'RNA', metric, count, 'read_pairs'])
    for name, fields, rows in [
        (f'{sample}.rna_export_retention.tsv', ['sample_id', 'modality', 'group', 'branch', 'metric', 'pairs', 'unit'], metrics),
        (f'{sample}.rna_export_pools.tsv', ['sample_id', 'group', 'split_id', 'pool_id', 'CB', 'XI', 'SB'], [pools[key] for key in sorted(pools)]),
    ]:
        with (output_dir / name).open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(fields)
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--prefix-length', type=int, required=True)
    parser.add_argument('--umi-length', type=int, required=True)
    parser.add_argument('--minimum-length', type=int, required=True)
    parser.add_argument('--pigz', type=Path, required=True)
    parser.add_argument('--threads', type=int, required=True, help='Total task CPU budget, shared by conversion and both mate compressors')
    args = parser.parse_args()
    if min(args.prefix_length, args.umi_length, args.minimum_length, args.threads) < 1:
        parser.error('Lengths and threads must be positive')
    export_fastqs(json.loads(args.manifest.read_text()), Path('.'), args.sample, args.prefix_length, args.umi_length, args.minimum_length, args.pigz, args.threads)


if __name__ == '__main__':
    main()
