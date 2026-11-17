#!/usr/bin/env python3
"""Compare `osmotool annotate` system calls on the reference genomes with
expected_systems.tsv, for one or more releases.

  python check_reference_systems.py --release releases/v9 [--release releases/v10 ...]

Reads <release>/real_genome_validation/reference_genomes/<genome>/<genome>.systems.tsv
(run_reference_genomes.sh) and prints, per expectation, the call in each
release, then every remaining system call and the per-genome family counts.
Writes reference_systems_check.tsv next to the first release's output.
Exit status 1 if any expectation fails in the LAST release given.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load_systems(release):
    base = os.path.join(release, 'real_genome_validation', 'reference_genomes')
    out = {}
    for genome in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        path = os.path.join(base, genome, f'{genome}.systems.tsv')
        if os.path.exists(path):
            out[genome] = {r['system']: r for r in csv.DictReader(open(path), delimiter='\t')}
    return out


def load_counts(release):
    base = os.path.join(release, 'real_genome_validation', 'reference_genomes')
    out = {}
    for genome in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        path = os.path.join(base, genome, f'{genome}.gene_counts.tsv')
        if not os.path.exists(path):
            continue
        counts = {}
        for line in open(path):
            if line.startswith('#') or line.startswith('gene\t'):
                continue
            f = line.rstrip('\n').split('\t')
            if float(f[1]) > 0:
                counts[f[0]] = int(float(f[1]))
        out[genome] = counts
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--release', action='append', required=True)
    ap.add_argument('--expected', default=os.path.join(HERE, 'expected_systems.tsv'))
    args = ap.parse_args()

    expected = [r for r in csv.DictReader(
        (l for l in open(args.expected) if not l.startswith('#')), delimiter='\t')]
    names = [os.path.basename(os.path.normpath(r)) for r in args.release]
    systems = [load_systems(r) for r in args.release]
    for name, s in zip(names, systems):
        if not s:
            sys.exit(f'no systems.tsv under {name}/real_genome_validation/reference_genomes/')

    def call(s, genome, system):
        row = s.get(genome, {}).get(system)
        if row is None:
            return None, 'missing'
        present = row['present'] == 'True'
        detail = f"{row['subunit_count']} subunits"
        if row['rescued_hits'] not in ('-', ''):
            detail += f", rescued {row['rescued_hits']}"
        return present, detail

    rows, n_fail = [], [0] * len(names)
    print(f'{"genome":<16}{"system":<16}{"expected":<10}' + ''.join(f'{n:<34}' for n in names))
    for e in expected:
        line = f'{e["genome"]:<16}{e["system"]:<16}{e["expected"]:<10}'
        rec = dict(genome=e['genome'], system=e['system'], expected=e['expected'], basis=e['basis'])
        for i, (name, s) in enumerate(zip(names, systems)):
            present, detail = call(s, e['genome'], e['system'])
            if e['expected'] == 'info':
                verdict = 'info'
            else:
                ok = present is not None and present == (e['expected'] == 'present')
                n_fail[i] += not ok
                verdict = 'ok' if ok else 'FAIL'
            cell = f'{verdict} ({"present" if present else "absent"}; {detail})'
            line += f'{cell:<34}'
            rec[name] = 'present' if present else 'absent'
            rec[name + '_ok'] = verdict
            rec[name + '_detail'] = detail
        print(line)
        rows.append(rec)
    print()
    n_scored = sum(1 for e in expected if e['expected'] != 'info')
    for name, nf in zip(names, n_fail):
        print(f'{name}: {n_scored - nf}/{n_scored} expectations met')
    for e in expected:
        if e['expected'] == 'info':
            print(f'  info  {e["genome"]} {e["system"]}: {e["note"]}')

    covered = {(e['genome'], e['system']) for e in expected}
    print('\nsystem calls with no stated expectation (present in any release):')
    for genome in sorted(set().union(*[set(s) for s in systems])):
        for system in sorted(set().union(*[set(s.get(genome, {})) for s in systems])):
            if (genome, system) in covered:
                continue
            calls = [call(s, genome, system) for s in systems]
            if any(c[0] for c in calls):
                print(f'  {genome:<16}{system:<16}' + ''.join(
                    f'{("present" if c[0] else "absent") + "; " + c[1]:<34}' for c in calls))

    counts = [load_counts(r) for r in args.release]
    print('\nfamilies detected per genome (n families; changes between first and last release):')
    for genome in sorted(counts[0]):
        a, b = counts[0][genome], counts[-1].get(genome, {})
        gained = sorted(set(b) - set(a)); lost = sorted(set(a) - set(b))
        changed = sorted(f for f in set(a) & set(b) if a[f] != b[f])
        line = f'  {genome:<16}' + '  '.join(f'{n}: {len(c.get(genome, {}))}' for n, c in zip(names, counts))
        if len(names) > 1:
            line += f'   gained {gained or "-"}  lost {lost or "-"}  copy number changed {changed or "-"}'
        print(line)

    out = os.path.join(args.release[-1], 'real_genome_validation', 'reference_systems_check.tsv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f'\nwrote {out}')
    sys.exit(1 if n_fail[-1] else 0)


if __name__ == '__main__':
    main()
