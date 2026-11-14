#!/usr/bin/env python3
"""Build the strain-level phenotype table from harvested BacDive records.

One row per strain that has BOTH a NaCl halophily record and a genome
accession. nacl_max_pct is the highest NaCl % with growth recorded positive
(tested relation 'growth' or 'maximum'); nogrowth_at_pct is the lowest % with
growth 'no' -- an upper bound on true tolerance, which makes the endpoint
interval-censored rather than a point value.

Also stamps each strain's osmoadaptation strategy from salt_in_taxa.tsv. That
assignment is pre-registered and taxonomy-based; see README.md.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re

LEVEL = {'complete': 0, 'chromosome': 1, 'scaffold': 2, 'contig': 3}
HERE = os.path.dirname(os.path.abspath(__file__))


def parse_conc(c: str):
    """'3-25 %' -> (3.0, 25.0, False); '<7 %' -> (7.0, 7.0, True); else None."""
    s = re.sub(r'\s+', '', str(c).replace('(w/v)', '').replace('%', ''))
    censored = s[:1] in ('<', '>')
    s = s.lstrip('<>')
    m = re.fullmatch(r'(\d*\.?\d+)-(\d*\.?\d+)', s)
    if m:
        return float(m.group(1)), float(m.group(2)), censored
    m = re.fullmatch(r'(\d*\.?\d+)', s)
    if m:
        return float(m.group(1)), float(m.group(1)), censored
    return None


def load_strategies(path):
    strat = {}
    if not os.path.exists(path):
        return strat
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('genus\t'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 4:
                strat[f[0]] = f[3]
    return strat


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True, help='output dir holding bacdive_raw/')
    ap.add_argument('--salt-in', default=os.path.join(HERE, 'salt_in_taxa.tsv'))
    args = ap.parse_args()

    strat = load_strategies(args.salt_in)
    rows, unparsed, molar = [], [], 0

    for f in sorted(glob.glob(os.path.join(args.dir, 'bacdive_raw', 'b*.json'))):
        # a batch with no strains returns results as [] rather than {}
        for bid, r in (json.load(open(f))['results'] or {}).items():
            h = (r.get('Physiology and metabolism') or {}).get('halophily')
            gs = (r.get('Sequence information') or {}).get('Genome sequences')
            if not h or not gs:
                continue
            gl = [g for g in (gs if isinstance(gs, list) else [gs])
                  if isinstance(g, dict) and g.get('INSDC accession')]
            if not gl:
                continue
            # prefer the most complete assembly, then the best strain match
            gl.sort(key=lambda g: (LEVEL.get(str(g.get('assembly level', '')).lower(), 9),
                                   -float(g.get('score') or 0)))
            g = gl[0]
            nt = r.get('Name and taxonomic classification') or {}
            mx = neg = None
            censored = False
            for e in (h if isinstance(h, list) else [h]):
                if not isinstance(e, dict) or e.get('salt') != 'NaCl':
                    continue
                c = str(e.get('concentration', ''))
                rel, gr = e.get('tested relation'), e.get('growth')
                if 'M' in c and '%' not in c:      # molar units, needs conversion
                    molar += 1
                    continue
                p = parse_conc(c)
                if p is None:
                    if c:
                        unparsed.append((bid, c))
                    continue
                if gr == 'positive' and rel in ('growth', 'maximum') and (mx is None or p[1] > mx):
                    mx, censored = p[1], p[2]
                if gr == 'no':
                    neg = p[0] if neg is None else min(neg, p[0])
            if mx is None and neg is None:
                continue
            asm = re.match(r'(\S+) assembly', str(g.get('description', '')))
            genus = nt.get('genus')
            rows.append(dict(
                bacdive_id=bid, species=nt.get('species'), genus=genus,
                family=nt.get('family'), domain=nt.get('domain'),
                strategy=strat.get(genus, 'salt_out'),
                nacl_max_pct=mx, censored=int(censored), nogrowth_at_pct=neg,
                insdc=g['INSDC accession'], asm_name=asm.group(1) if asm else None,
                asm_level=g.get('assembly level'), match_score=g.get('score')))

    if not rows:
        raise SystemExit(f'no usable strains found under {args.dir}/bacdive_raw/')

    rows.sort(key=lambda r: (r['nacl_max_pct'] is None, r['nacl_max_pct'] or 0))
    out = os.path.join(args.dir, 'pilot_strains.tsv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    n_quant = sum(1 for r in rows if r['nacl_max_pct'] is not None)
    genera = {r['genus'] for r in rows}
    multi = sum(1 for g in genera
                if sum(1 for r in rows if r['genus'] == g) >= 3)
    print(f'strains with phenotype + genome : {len(rows)}')
    print(f'  quantitative NaCl max         : {n_quant}')
    print(f'  no-growth upper bound         : {sum(1 for r in rows if r["nogrowth_at_pct"] is not None)}')
    print(f'  salt_in (pre-registered)      : {sum(1 for r in rows if r["strategy"] == "salt_in")}')
    print(f'  archaea                       : {sum(1 for r in rows if r["domain"] == "Archaea")}')
    print(f'  distinct genera               : {len(genera)}')
    print(f'  genera with >=3 strains (SetB): {multi}')
    print(f'  molar-unit records skipped    : {molar}')
    print(f'  unparsed concentrations       : {len(unparsed)} {unparsed[:5]}')
    print(f'wrote                           : {out}')


if __name__ == '__main__':
    main()
