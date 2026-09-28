#!/usr/bin/env python3
"""Small end-to-end barcode prior test; no production data or timing assertions."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chromap', required=True, type=Path)
    parser.add_argument('--contract-runner', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    root = args.out.resolve()
    commands = []

    def run(name, argv):
        commands.append(dict(name=name, argv=[str(x) for x in argv]))
        (root / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        result = subprocess.run(argv, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (root / (name + '.log')).write_text(result.stdout)
        result.check_returncode()
        return result.stdout

    rng = random.Random(17)
    genome = ''.join(rng.choice('ACGT') for _ in range(50000))
    (root / 'ref.fa').write_text('>chrTest\n' + genome + '\n')
    a, b, ambiguous = 'ACGT' * 4, 'ACGT' * 3 + 'ACGA', 'ACGT' * 3 + 'ACGC'
    (root / 'whitelist.txt').write_text(a + '\n' + b + '\n')
    # Early priors favor A; complete-input priors favor B. The ambiguous
    # barcode is one mismatch from both and clears 0.9 under either model.
    lanes = [[a] * 12 + [b], [b] * 121 + [ambiguous]]
    ordinal = 0
    for lane, barcodes in enumerate(lanes):
        streams = {kind: [] for kind in ('r1', 'r2', 'bc')}
        for barcode in barcodes:
            start = 200 + ordinal * 250
            seqs = dict(r1=genome[start:start + 90],
                        r2=genome[start + 100:start + 190].translate(str.maketrans('ACGT', 'TGCA'))[::-1],
                        bc=barcode)
            for kind, seq in seqs.items():
                streams[kind].append(f'@read{ordinal}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
            ordinal += 1
        for kind, records in streams.items():
            (root / f'{kind}.{lane}.fq').write_text(''.join(records))
    run('index', [str(args.chromap), '--build-index', '-r', str(root / 'ref.fa'),
                  '-o', str(root / 'ref.idx'), '-k', '11', '-w', '5'])
    csv = lambda kind: ','.join(str(root / f'{kind}.{lane}.fq') for lane in range(2))
    routes = [('single', args.chromap, False), ('paired', args.chromap, True)]
    if args.contract_runner:
        routes.append(('contract', args.contract_runner, True))
    results = []
    for route, binary, paired in routes:
        outputs = {}
        for mode, limit in [('default', None), ('explicit', 20000000), ('bounded', 1), ('all', 0)]:
            name = route + '_' + mode
            output = root / (name + '.bed')
            argv = [str(binary), '--ref', str(root / 'ref.fa'), '--index', str(root / 'ref.idx'),
                    '--read1', csv('r1'), '--barcode', csv('bc'),
                    '--barcode-whitelist', str(root / 'whitelist.txt'), '--output', str(output)]
            if paired:
                argv += ['--read2', csv('r2')]
            if limit is not None:
                argv += ['--barcode-sample-limit', str(limit)]
            log = run(name, argv)
            learned = int(re.search(r'Compute barcode abundance using (\d+)', log).group(1))
            assert learned == (13 if mode == 'bounded' else 134), (name, learned)
            if mode in ('default', 'explicit'):
                assert 'limit=20000000 exact whitelist observations' in log
            if mode == 'all':
                assert 'Barcode abundance sampling: all barcode inputs.' in log
            counts = Counter(line.split('\t')[3] for line in output.read_text().splitlines()
                             if line and not line.startswith('#'))
            expected = Counter({a: 13, b: 122}) if mode == 'bounded' else Counter({a: 12, b: 123})
            assert counts == expected, (name, counts, expected)
            outputs[mode] = output.read_bytes()
            results.append(dict(route=route, mode=mode, learned=learned, mapped_fragments=sum(counts.values())))
        assert outputs['default'] == outputs['explicit'] == outputs['all'], route
    binaries = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [args.chromap, args.contract_runner] if p is not None}
    (root / 'results.json').write_text(json.dumps(dict(status='passed', binaries=binaries, cases=results), indent=2) + '\n')
    print(f'PASS: {len(results)} barcode sampling cases; all 135 reads mapped; ambiguous correction follows learned priors')


if __name__ == '__main__':
    main()
