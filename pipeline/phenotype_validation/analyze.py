#!/usr/bin/env python3
"""Read out the phenotype validation against the prespecified plan in README.md
(see its "Plan provenance" section: a plan committed to git before the full
data were collected, not an externally registered protocol).

Structure follows that plan rather than exploring:

  Set A   unbiased estimate, panel always compared against its baselines, with
          BOTH a strain bootstrap and a genus-cluster bootstrap. The cluster
          version is the honest one -- strains of a genus are not independent.
  Set B   within-genus contrast, the only part that controls for ancestry.
  salt-in prespecified negative control, with a family-specific prediction
          (K+ uptake present, compatible solutes absent) -- NOT "low breadth".
  power   from a PRE-SPECIFIED effect size. Never from an observed correlation:
          that is circular, and biased upward because a small sample that looks
          impressive is usually an overestimate.

No scipy: Spearman, partial correlation and both bootstraps are implemented
here so this runs outside the container too.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import statistics
import sys

# Prespecified family groups for the salt-in prediction (README.md).
# Salt-in organisms accumulate molar KCl, so bulk K+ uptake IS their mechanism
# and should be present; the compatible-solute machinery should not be.
K_UPTAKE = ['trkA', 'trkH', 'ktrA', 'ktrB', 'ktrD']
SOLUTE = ['ectA', 'ectB', 'ectC', 'betL', 'betA', 'betB', 'otsA', 'otsB',
          'proX', 'proP', 'opuAA', 'opuAB', 'opuAC', 'opuBA', 'opuBB', 'opuBC',
          'opuCA', 'opuCB', 'opuCC']

# Smallest correlation worth acting on -> required n. Chosen independently of
# any observed result; see README.md "Power".
PRESPECIFIED_RHO = [0.30, 0.20, 0.15, 0.10]

META = {'bacdive_id', 'species', 'genus', 'family', 'domain', 'strategy',
        'nacl_max_pct', 'censored', 'nogrowth_at_pct', 'asm_level',
        'total_proteins', 'n_families_present', 'total_osmo_copies'}


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float('nan')
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    if sx == 0 or sy == 0:
        return float('nan')
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def partial_spearman(x, y, z):
    """Rank-level correlation of x,y with z regressed out of both."""
    rx, ry, rz = ranks(x), ranks(y), ranks(z)

    def resid(v):
        b = pearson(rz, v)
        if math.isnan(b):
            return v
        mz, mv = sum(rz) / len(rz), sum(v) / len(v)
        sz = statistics.pstdev(rz) or 1.0
        sv = statistics.pstdev(v) or 1.0
        slope = b * sv / sz
        return [vi - (mv + slope * (zi - mz)) for vi, zi in zip(v, rz)]

    return pearson(resid(rx), resid(ry))


def ci(vals):
    vals = sorted(v for v in vals if not math.isnan(v))
    if len(vals) < 100:
        return float('nan'), float('nan')
    return vals[int(.025 * len(vals))], vals[int(.975 * len(vals))]


def boot_strain(fn, cols, n=4000, rng=None):
    rng = rng or random
    N = len(cols[0])
    out = []
    for _ in range(n):
        idx = [rng.randrange(N) for _ in range(N)]
        try:
            out.append(fn(*[[c[i] for i in idx] for c in cols]))
        except Exception:
            pass
    return ci(out)


def boot_cluster(fn, cols, groups, n=4000, rng=None):
    """Resample whole genera, not strains: strains of a genus aren't independent."""
    rng = rng or random
    by = {}
    for i, g in enumerate(groups):
        by.setdefault(g, []).append(i)
    keys = list(by)
    out = []
    for _ in range(n):
        idx = []
        for _ in range(len(keys)):
            idx.extend(by[keys[rng.randrange(len(keys))]])
        try:
            out.append(fn(*[[c[i] for i in idx] for c in cols]))
        except Exception:
            pass
    return ci(out)


def n_for_rho(rho):
    """Sample size for 80% power at alpha .05 to detect a correlation of rho."""
    z = abs(math.atanh(abs(rho)))
    return (1.959964 + 0.841621) ** 2 / (z * z) + 3


def eta_sq(values, groups):
    """Fraction of variance in `values` explained by `groups` alone."""
    gm = statistics.mean(values)
    total = sum((v - gm) ** 2 for v in values)
    if total == 0:
        return float('nan')
    by = {}
    for v, g in zip(values, groups):
        by.setdefault(g, []).append(v)
    between = sum(len(vs) * (statistics.mean(vs) - gm) ** 2 for vs in by.values())
    return between / total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--min-genus-n', type=int, default=3, help='Set B: strains per genus')
    ap.add_argument('--min-genus-spread', type=float, default=4.0,
                    help='Set B: required NaCl %% spread within a genus')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    path = os.path.join(args.dir, 'feature_matrix.tsv')
    if not os.path.exists(path):
        sys.exit(f'missing {path} -- run build_matrix.py first')
    allrows = [r for r in csv.DictReader(open(path), delimiter='\t') if r['nacl_max_pct']]
    fams = [c for c in allrows[0] if c not in META]

    salt_in = [r for r in allrows if r.get('strategy') == 'salt_in']
    rows = [r for r in allrows if r.get('strategy') != 'salt_in']

    def col(rs, name):
        return [float(r[name]) for r in rs]

    y = col(rows, 'nacl_max_pct')
    prot = col(rows, 'total_proteins')
    breadth = col(rows, 'n_families_present')
    copies = col(rows, 'total_osmo_copies')
    genera = [r['genus'] for r in rows]

    print('=' * 72)
    print('SET A -- primary, salt-out strains only (prespecified)')
    print('=' * 72)
    print(f'n strains : {len(rows)}   n genera : {len(set(genera))}')
    print(f'NaCl max  : {min(y):.1f}-{max(y):.1f} % (median {statistics.median(y):.1f})')
    print(f'salt-in strains held out for the negative control : {len(salt_in)}')
    print()
    print('Spearman rho vs NaCl max   [95% strain boot] [95% GENUS-CLUSTER boot]')
    for label, v in (('panel breadth', breadth),
                     ('panel total copies', copies),
                     ('BASELINE total_proteins', prot)):
        rho = spearman(v, y)
        a, b = boot_strain(spearman, [v, y], rng=rng)
        c, d = boot_cluster(spearman, [v, y], genera, rng=rng)
        print(f'  {label:<24} {rho:+.3f}  [{a:+.3f},{b:+.3f}]  [{c:+.3f},{d:+.3f}]')

    print()
    print('with genome size partialled out')
    for label, v in (('breadth | total_proteins', breadth),
                     ('copies | total_proteins', copies)):
        rho = partial_spearman(v, y, prot)
        a, b = boot_strain(partial_spearman, [v, y, prot], rng=rng)
        c, d = boot_cluster(partial_spearman, [v, y, prot], genera, rng=rng)
        print(f'  {label:<24} {rho:+.3f}  [{a:+.3f},{b:+.3f}]  [{c:+.3f},{d:+.3f}]')

    print()
    print('confound, quantified: fraction of NaCl variance explained by taxonomy alone')
    per_genus = len(rows) / max(len(set(genera)), 1)
    for lvl in ('genus', 'family'):
        if lvl in rows[0]:
            n_grp = len(set(r[lvl] for r in rows))
            e = eta_sq(y, [r[lvl] for r in rows])
            flag = '  <-- DEGENERATE, ignore' if len(rows) / max(n_grp, 1) < 2 else ''
            print(f'  eta^2 by {lvl:<7} = {e:.3f}  ({n_grp} groups, '
                  f'{len(rows)/max(n_grp,1):.1f} strains/group){flag}')
    if per_genus < 2:
        print()
        print(f'  WARNING: {per_genus:.1f} strains per genus. With ~1 strain per group,')
        print('  between-group variance IS total variance, so eta^2 is forced toward')
        print('  1.0 by construction and says nothing about phylogenetic signal. The')
        print('  genus-cluster bootstrap is likewise degenerate -- resampling clusters')
        print('  of size 1 reproduces the strain bootstrap, which is why the two CIs')
        print('  above agree. All THREE ancestry controls (cluster boot, eta^2, Set B)')
        print('  only become meaningful after the targeted multi-strain-per-genus')
        print('  harvest. Treat Set A as uncontrolled for phylogeny until then.')

    print()
    print('=' * 72)
    print('SET B -- within-genus contrast (the only ancestry-controlled test)')
    print('=' * 72)
    by_genus = {}
    for r in rows:
        by_genus.setdefault(r['genus'], []).append(r)
    qualify = {g: rs for g, rs in by_genus.items()
               if len(rs) >= args.min_genus_n
               and (max(col(rs, 'nacl_max_pct')) - min(col(rs, 'nacl_max_pct')))
               >= args.min_genus_spread}
    print(f'genera with >={args.min_genus_n} strains spanning >={args.min_genus_spread} pp NaCl : {len(qualify)}')
    if not qualify:
        print()
        print('  SET B IS EMPTY -- no ancestry-controlled estimate is possible from')
        print('  this sample. A random BacDive draw yields roughly one strain per')
        print('  genus, so Set B requires the targeted harvest described in README.md.')
        print('  Any Set A result below remains fully confounded with phylogeny.')
    else:
        within = []
        for g, rs in sorted(qualify.items()):
            rho = spearman(col(rs, 'n_families_present'), col(rs, 'nacl_max_pct'))
            within.append(rho)
            print(f'  {g:<22} n={len(rs):<3} rho={rho:+.3f}')
        good = [w for w in within if not math.isnan(w)]
        if good:
            print(f'  median within-genus rho = {statistics.median(good):+.3f} '
                  f'over {len(good)} genera')
            print('  DIAGNOSTIC: positive here => signal survives ancestry control.')
            print('              collapsing to ~0 while Set A stays high => it was ancestry.')

    print()
    print('=' * 72)
    print('SALT-IN NEGATIVE CONTROL (prespecified, family-specific)')
    print('=' * 72)
    if not salt_in:
        print('  no salt-in strains present -- needs the targeted /v2/taxon harvest')
        print('  over the genera in salt_in_taxa.tsv. Until then the panel has not')
        print('  been tested against a high-phenotype/mechanism-absent cohort.')
    else:
        def grp_mean(rs, group):
            present = [c for c in group if c in fams]
            return statistics.mean(
                sum(1 for c in present if float(r[c]) > 0) / len(present) for r in rs)
        print(f'  salt-in n={len(salt_in)}, salt-out n={len(rows)}')
        print(f'  {"group":<28}{"salt-in":>10}{"salt-out":>10}  prediction')
        print(f'  {"K+ uptake (trk/ktr)":<28}{grp_mean(salt_in, K_UPTAKE):>10.2f}'
              f'{grp_mean(rows, K_UPTAKE):>10.2f}  salt-in present/high')
        print(f'  {"compatible solutes":<28}{grp_mean(salt_in, SOLUTE):>10.2f}'
              f'{grp_mean(rows, SOLUTE):>10.2f}  salt-in low/absent')
        print('  (fraction of each family group present, averaged over strains)')
        print('  Panel firing on compatible solutes in salt-in strains is a FALSE')
        print('  POSITIVE no other stratum would reveal.')

    print()
    print('=' * 72)
    print('POWER -- from pre-specified effect sizes, not from any result above')
    print('=' * 72)
    for rho in PRESPECIFIED_RHO:
        print(f'  to detect rho={rho:.2f} at 80% power, alpha .05 : n ~ {n_for_rho(rho):.0f}')
    print(f'  current usable n = {len(rows)}')


if __name__ == '__main__':
    main()
