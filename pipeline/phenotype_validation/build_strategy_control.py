#!/usr/bin/env python3
"""Build the strain table for the within-domain salt-in control.

The salt-in control in analyze.py compares salt-in strains with all salt-out
strains. In the NaCl-phenotyped set nearly every salt-in strain is an
archaeon and nearly every salt-out strain a bacterium, and the 29 salt-out
archaea available there score as low on compatible-solute families as the
salt-in strains do. That comparison cannot tell osmotic strategy from domain.

The gene-content prediction does not need a NaCl measurement, only a genome
and a lineage. This script therefore takes every BacDive strain with a genome
accession, whether or not it has a halophily record, from:

  archaea            all of them: class Halobacteria (salt-in) against every
                     other archaeal class (salt-out)
  salt-in bacteria   order Halanaerobiales, genera Salinibacter and Salisaeta

and writes <dir>/strategy_control/pilot_strains.tsv in build_table.py's
format (plus <dir>/bacdive_classes.tsv, the class of every BacDive strain), so fetch_genomes.sh, run_random_panel.sh and run_annotate_cached.sh
run on that directory unchanged. Genomes already downloaded for the main
analysis are linked, not fetched again. The salt-out bacterial reference is
the main feature matrix; analyze_strategy_control.py joins the two.

Strategy is assigned from lineage, by the criteria in salt_in_taxa.tsv's
header. The genus list in that file is checked against the lineage rule and
any disagreement is printed.

This comparison was added after the main results were known (2026-11-14).
It is a post hoc extension of the control, not part of the original plan.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re

from build_table import LEVEL, load_strategies, parse_conc

HERE = os.path.dirname(os.path.abspath(__file__))
SALT_IN_CLASSES = {'Halobacteria'}
SALT_IN_ORDERS = {'Halanaerobiales'}
SALT_IN_FAMILIES = {'Halanaerobiaceae', 'Halobacteroidaceae'}   # order sometimes unassigned
SALT_IN_GENERA = {'Salinibacter', 'Salisaeta'}


def lineage_strategy(nt: dict) -> str:
    if (nt.get('class') in SALT_IN_CLASSES or nt.get('order') in SALT_IN_ORDERS
            or nt.get('family') in SALT_IN_FAMILIES or nt.get('genus') in SALT_IN_GENERA):
        return 'salt_in'
    return 'salt_out'


def nacl_max(h):
    """Highest NaCl % with positive growth, as in build_table.py; None if absent."""
    mx = None
    for e in (h if isinstance(h, list) else [h]) if h else []:
        if not isinstance(e, dict) or e.get('salt') != 'NaCl':
            continue
        c = str(e.get('concentration', ''))
        if 'M' in c and '%' not in c:
            continue
        p = parse_conc(c)
        if p and e.get('growth') == 'positive' and e.get('tested relation') in ('growth', 'maximum'):
            mx = p[1] if mx is None else max(mx, p[1])
    return mx


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True, help='main output dir holding bacdive_raw/ and genomes/')
    ap.add_argument('--salt-in', default=os.path.join(HERE, 'salt_in_taxa.tsv'))
    args = ap.parse_args()

    listed = load_strategies(args.salt_in)
    out_dir = os.path.join(args.dir, 'strategy_control')
    os.makedirs(os.path.join(out_dir, 'genomes'), exist_ok=True)

    rows, disagree, classes = [], [], []
    for f in sorted(glob.glob(os.path.join(args.dir, 'bacdive_raw', 'b*.json'))):
        for bid, r in (json.load(open(f))['results'] or {}).items():
            nt = r.get('Name and taxonomic classification') or {}
            classes.append((bid, nt.get('class') or ''))
            strategy = lineage_strategy(nt)
            by_list = listed.get(nt.get('genus'), 'salt_out')
            if nt.get('domain') != 'Archaea' and strategy != 'salt_in':
                if by_list == 'salt_in' and nt.get('domain') == 'Bacteria':
                    disagree.append((bid, nt.get('genus'), 'listed salt_in, lineage rule salt_out'))
                continue
            if by_list != strategy:
                disagree.append((bid, nt.get('genus'), f'list {by_list}, lineage rule {strategy}'))
            gs = (r.get('Sequence information') or {}).get('Genome sequences')
            gl = [g for g in (gs if isinstance(gs, list) else [gs] if gs else [])
                  if isinstance(g, dict) and g.get('INSDC accession')]
            if not gl:
                continue
            gl.sort(key=lambda g: (LEVEL.get(str(g.get('assembly level', '')).lower(), 9),
                                   -float(g.get('score') or 0)))
            g = gl[0]
            asm = re.match(r'(\S+) assembly', str(g.get('description', '')))
            mx = nacl_max((r.get('Physiology and metabolism') or {}).get('halophily'))
            rows.append(dict(
                bacdive_id=bid, species=nt.get('species'), genus=nt.get('genus'),
                family=nt.get('family'), order=nt.get('order'), tax_class=nt.get('class'),
                domain=nt.get('domain'), strategy=strategy,
                nacl_max_pct='' if mx is None else mx,
                insdc=g['INSDC accession'], asm_name=asm.group(1) if asm else None,
                asm_level=g.get('assembly level'), match_score=g.get('score')))

    rows.sort(key=lambda r: (r['domain'] or '', r['strategy'], r['tax_class'] or '', int(r['bacdive_id'])))
    out = os.path.join(out_dir, 'pilot_strains.tsv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    # class of every BacDive strain, for class-matched references
    with open(os.path.join(args.dir, 'bacdive_classes.tsv'), 'w') as fh:
        fh.write('bacdive_id\ttax_class\n')
        for bid, cls in classes:
            fh.write(f'{bid}\t{cls}\n')

    # link genomes the main analysis already holds
    have = {}
    main_g = os.path.join(args.dir, 'genomes')
    for name in os.listdir(main_g) if os.path.isdir(main_g) else []:
        if name.endswith(('.fna', '.fna.gz')):
            have[name.split('_', 1)[0]] = name
    linked = 0
    for r in rows:
        name = have.get(r['bacdive_id'])
        if name:
            dst = os.path.join(out_dir, 'genomes', name)
            if not os.path.lexists(dst):
                os.symlink(os.path.join('..', '..', 'genomes', name), dst)
            linked += 1

    def n(pred):
        return sum(1 for r in rows if pred(r))
    print(f'strains with a genome            : {len(rows)}')
    print(f'  archaea, salt-in (Halobacteria): {n(lambda r: r["domain"] == "Archaea" and r["strategy"] == "salt_in")}')
    print(f'  archaea, salt-out              : {n(lambda r: r["domain"] == "Archaea" and r["strategy"] == "salt_out")}')
    print(f'  bacteria, salt-in              : {n(lambda r: r["domain"] == "Bacteria" and r["strategy"] == "salt_in")}')
    print(f'  with a quantitative NaCl max   : {n(lambda r: r["nacl_max_pct"] != "")}')
    print(f'genomes linked from main analysis: {linked}  (to fetch: {len(rows) - linked})')
    print(f'genus list vs lineage rule, strains that disagree: {len(disagree)}')
    for bid, genus, what in disagree[:60]:
        print(f'  {bid:>8}  {genus}: {what}')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
