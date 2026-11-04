#!/usr/bin/env python3
"""Join osmotool annotate output to the BacDive phenotype table.

Reads <dir>/annotate/<bid>/<bid>.gene_counts.tsv for every annotated strain and
writes <dir>/feature_matrix.tsv: phenotype columns from pilot_strains.tsv, then
per-family raw counts, plus the breadth/baseline summaries the analysis
compares against each other.

No pandas: the matrix is small and the scientific stack is not guaranteed to be
present outside the container.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

PHENO_COLS = ['bacdive_id', 'species', 'genus', 'family', 'domain', 'strategy',
              'nacl_max_pct', 'censored', 'nogrowth_at_pct', 'asm_level']


def read_counts(path):
    """-> (total_proteins, {family: raw_count})"""
    total, counts = None, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('# total_proteins'):
                total = int(line.split('\t')[1])
            elif line.startswith('#') or line.startswith('gene\t'):
                continue
            else:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 2:
                    counts[f[0]] = float(f[1])
    return total, counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    args = ap.parse_args()

    tsv = os.path.join(args.dir, 'pilot_strains.tsv')
    if not os.path.exists(tsv):
        sys.exit(f'missing {tsv} -- run build_table.py first')
    pheno = {r['bacdive_id']: r for r in csv.DictReader(open(tsv), delimiter='\t')}

    rows, families = [], set()
    for path in sorted(glob.glob(os.path.join(args.dir, 'annotate', '*', '*.gene_counts.tsv'))):
        bid = os.path.basename(path).split('.')[0]
        if bid not in pheno:
            continue
        total, counts = read_counts(path)
        families |= set(counts)
        rows.append((bid, total, counts))

    if not rows:
        sys.exit(f'no gene_counts.tsv under {args.dir}/annotate/ -- run run_annotate.sh first')

    families = sorted(families)
    out = os.path.join(args.dir, 'feature_matrix.tsv')
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(PHENO_COLS + ['total_proteins', 'n_families_present',
                                 'total_osmo_copies'] + families)
        for bid, total, counts in rows:
            p = pheno[bid]
            present = sum(1 for f in families if counts.get(f, 0) > 0)
            w.writerow([p.get(c, '') for c in PHENO_COLS]
                       + [total, present, f'{sum(counts.get(f, 0) for f in families):.1f}']
                       + [f'{counts.get(f, 0):.1f}' for f in families])

    breadths = sorted(sum(1 for f in families if c.get(f, 0) > 0) for _, _, c in rows)
    print(f'strains annotated      : {len(rows)}')
    print(f'families in output     : {len(families)}')
    print(f'families present/strain: min {breadths[0]} median '
          f'{breadths[len(breadths)//2]} max {breadths[-1]}')
    print(f'wrote                  : {out}')


if __name__ == '__main__':
    main()
