#!/usr/bin/env python3
"""Join osmotool annotate output to the BacDive phenotype table.

Reads <dir>/annotate/<bid>/<bid>.gene_counts.tsv for every annotated strain and
writes <dir>/feature_matrix.tsv: phenotype columns from pilot_strains.tsv, then
per-family raw counts, plus the breadth/baseline summaries the analysis
compares against each other.

Assemblies under --min-genome-bp (default 1 Mb of A/C/G/T) are left out of the
matrix and listed in the log. BacDive links a few strains to partial
assemblies (a single 5.3 kb scaffold, a 68 kb set of contigs), which carry
almost no panel families and would pass for low-breadth genomes. The smallest
complete prokaryote genomes in the sweep are ~1.3 Mb.

No pandas: the matrix is small and the scientific stack is not guaranteed to be
present outside the container.
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
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


# Only assemblies that could plausibly fall under the cutoff get their bases
# counted. The test uses the uncompressed size: a plain file's size on disk, or
# a gzip file's ISIZE trailer (its last 4 bytes, the uncompressed length of the
# final member). FASTA overhead is ~1%, so 2x the cutoff leaves a wide margin.
# A multi-member gzip under-reports in ISIZE, which only means its bases get
# counted exactly.
def uncompressed_size(path):
    if not path.endswith('.gz'):
        return os.path.getsize(path)
    with open(path, 'rb') as fh:
        fh.seek(-4, os.SEEK_END)
        return int.from_bytes(fh.read(4), 'little')


def genome_files(gdir):
    """{bacdive_id: path} for <gdir>/<bid>_<accession>.fna[.gz]"""
    out = {}
    for name in os.listdir(gdir) if os.path.isdir(gdir) else []:
        if name.endswith(('.fna', '.fna.gz')):
            out.setdefault(name.split('_', 1)[0], os.path.join(gdir, name))
    return out


def genome_bp(path, min_bp):
    """A/C/G/T count of a FASTA (.gz or plain). Files too large to fall under
    min_bp are returned as their byte size, uncounted."""
    size = uncompressed_size(path)
    if size >= 2 * min_bp:
        return size
    op = gzip.open if path.endswith('.gz') else open
    n = 0
    with op(path, 'rt') as fh:
        for line in fh:
            if not line.startswith('>'):
                s = line.upper()
                n += s.count('A') + s.count('C') + s.count('G') + s.count('T')
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', default=None,
                    help='a phenotype_validation directory; or give its parts with the options below')
    ap.add_argument('--genomes', default=None, help='default: <dir>/genomes')
    ap.add_argument('--pilot', default=None, help='default: <dir>/pilot_strains.tsv')
    ap.add_argument('--annotate', default=None, help='per-strain annotate output, <bid>/; default: <dir>/annotate')
    ap.add_argument('--out', default=None, help='default: <dir>/feature_matrix.tsv')
    ap.add_argument('--min-genome-bp', type=int, default=1_000_000)
    args = ap.parse_args()
    if args.dir is None and not (args.genomes and args.pilot and args.annotate and args.out):
        ap.error('give --dir, or --genomes, --pilot, --annotate and --out')
    gdir = args.genomes or os.path.join(args.dir, 'genomes')
    gfiles = genome_files(gdir)

    tsv = args.pilot or os.path.join(args.dir, 'pilot_strains.tsv')
    if not os.path.exists(tsv):
        sys.exit(f'missing {tsv} -- run build_table.py first')
    pheno = {r['bacdive_id']: r for r in csv.DictReader(open(tsv), delimiter='\t')}

    rows, families, small, nogenome = [], set(), [], 0
    annot = args.annotate or os.path.join(args.dir, 'annotate')
    for path in sorted(glob.glob(os.path.join(annot, '*', '*.gene_counts.tsv'))):
        bid = os.path.basename(path).split('.')[0]
        if bid not in pheno:
            continue
        if bid not in gfiles:
            nogenome += 1
        elif (bp := genome_bp(gfiles[bid], args.min_genome_bp)) < args.min_genome_bp:
            small.append((bid, bp))
            continue
        total, counts = read_counts(path)
        families |= set(counts)
        rows.append((bid, total, counts))

    if not rows:
        sys.exit(f'no gene_counts.tsv under {annot}/ -- run BACDIVE_ANNOTATE first')

    families = sorted(families)
    out = args.out or os.path.join(args.dir, 'feature_matrix.tsv')
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
    print(f'excluded, assembly < {args.min_genome_bp:,} bp: {len(small)}')
    for bid, bp in small:
        p = pheno[bid]
        print(f'  {bid:>8}  {bp:>9,} bp  {p.get("species", "")}  ({p.get("asm_level", "")})')
    if nogenome:
        print(f'WARNING: {nogenome} strains had no genome file in {gdir}; size not checked')
    print(f'strains in matrix      : {len(rows)}')
    print(f'families in output     : {len(families)}')
    print(f'families present/strain: min {breadths[0]} median '
          f'{breadths[len(breadths)//2]} max {breadths[-1]}')
    print(f'wrote                  : {out}')


if __name__ == '__main__':
    main()
