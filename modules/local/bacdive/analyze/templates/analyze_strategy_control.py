#!/usr/bin/env python3
"""Within-domain salt-in control (post hoc; see build_strategy_control.py).

Question: do salt-in organisms show the gene pattern their physiology
predicts -- K+ uptake families present, compatible-solute families scarce --
when they are compared with organisms of the SAME domain, so that the
contrast is not archaea against bacteria?

Comparisons:
  archaea    class Halobacteria (salt-in) vs all other archaea (salt-out)
  bacteria   Halanaerobiales + Salinibacter/Salisaeta (salt-in) vs
             (a) all salt-out bacteria of the main analysis,
             (b) those with NaCl max >= --halophile-pct (halophile-matched),
             (c) salt-out strains of the same taxonomic class

Per genome: families detected among the 5 K+ uptake families, among the 19
compatible-solute families, and among the 174 random pool families (how many
families of any kind the HMMs find in that genome). Strains of a genus are
not independent, so tests run on genus means (Mann-Whitney U) and intervals
come from a genus-cluster bootstrap. An OLS with genus-clustered standard
errors adjusts the group difference for the random-pool count and genome
size, i.e. for general detectability.

Usage:
  python analyze_strategy_control.py --dir <main PV dir> \\
      [--osmo-dir <dir with feature_matrix.tsv and annotate/ for this release>]
--dir holds strategy_control/ (strain table, random_panel/) and
random_panel/strains/ of the main analysis. --osmo-dir (default --dir) holds
the release-specific annotate output: <osmo-dir>/feature_matrix.tsv for the
main strains and <osmo-dir>/strategy_control/annotate/ for the control set.
Writes <osmo-dir>/strategy_control/{strategy_groups,strategy_family_prevalence}.tsv
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

K_UPTAKE = ['trkA', 'trkH', 'ktrA', 'ktrB', 'ktrD']
SOLUTE = ['ectA', 'ectB', 'ectC', 'betL', 'betA', 'betB', 'otsA', 'otsB',
          'proX', 'proP', 'opuAA', 'opuAB', 'opuAC', 'opuBA', 'opuBB', 'opuBC',
          'opuCA', 'opuCB', 'opuCC']
MIN_GENOME_BP = 1_000_000


def read_counts(path):
    counts = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('gene\t'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2:
                counts[f[0]] = float(f[1])
    return counts


def pool_stats(strain_dir, bid):
    """(n pool families detected, genome_bp, gc) or None if not computed."""
    pc = os.path.join(strain_dir, bid, f'{bid}.pool_counts.tsv')
    gs = os.path.join(strain_dir, bid, f'{bid}.genome_stats.tsv')
    if not (os.path.exists(pc) and os.path.exists(gs)):
        return None
    n_pool = max(sum(1 for _ in open(pc)) - 1, 0)
    g = pd.read_csv(gs, sep='\t').iloc[0]
    return n_pool, float(g['genome_bp']), float(g['gc'])


def load(args):
    sc = os.path.join(args.dir, 'strategy_control') if args.dir else None
    sc_osmo = os.path.join(args.osmo_dir, 'strategy_control') if args.osmo_dir else None
    ctrl_annot = args.control_annotate or os.path.join(sc_osmo, 'annotate')
    ctrl_strains = args.control_strains or os.path.join(sc, 'random_panel', 'strains')
    strains = args.strains or os.path.join(args.dir, 'random_panel', 'strains')
    rows = []

    ctrl = pd.read_csv(args.control_pilot or os.path.join(sc, 'pilot_strains.tsv'), sep='\t', dtype=str)
    missing = 0
    for r in ctrl.itertuples():
        gc_path = os.path.join(ctrl_annot, r.bacdive_id, f'{r.bacdive_id}.gene_counts.tsv')
        ps = pool_stats(ctrl_strains, r.bacdive_id)
        if not os.path.exists(gc_path) or ps is None:
            missing += 1
            continue
        rows.append(dict(bacdive_id=r.bacdive_id, genus=r.genus, tax_class=r.tax_class,
                         domain=r.domain, strategy=r.strategy,
                         nacl=float(r.nacl_max_pct) if isinstance(r.nacl_max_pct, str) and r.nacl_max_pct else np.nan,
                         n_pool=ps[0], genome_bp=ps[1], gc=ps[2], **read_counts(gc_path)))
    print(f'control set: {len(ctrl)} strains in table, {missing} without annotate/pool output')

    # salt-out bacteria: the main analysis matrix
    fm = pd.read_csv(args.matrix or os.path.join(args.osmo_dir, 'feature_matrix.tsv'), sep='\t',
                     dtype={'bacdive_id': str})
    fm = fm[(fm['domain'] == 'Bacteria') & (fm['strategy'] != 'salt_in')]
    classes = {}
    raw_tbl = args.classes or os.path.join(args.dir, 'bacdive_classes.tsv')
    if os.path.exists(raw_tbl):
        classes = dict(pd.read_csv(raw_tbl, sep='\t', dtype=str).values)
    meta_cols = {'bacdive_id', 'species', 'genus', 'family', 'domain', 'strategy', 'nacl_max_pct',
                 'censored', 'nogrowth_at_pct', 'asm_level', 'total_proteins',
                 'n_families_present', 'total_osmo_copies'}
    fams = [c for c in fm.columns if c not in meta_cols]
    have = {r['bacdive_id'] for r in rows}
    for r in fm.itertuples(index=False):
        d = r._asdict()
        if d['bacdive_id'] in have:
            continue
        ps = pool_stats(strains, d['bacdive_id'])
        if ps is None:
            continue
        rows.append(dict(bacdive_id=d['bacdive_id'], genus=d['genus'],
                         tax_class=classes.get(d['bacdive_id']), domain='Bacteria',
                         strategy='salt_out', nacl=d['nacl_max_pct'],
                         n_pool=ps[0], genome_bp=ps[1], gc=ps[2],
                         **{f: d[f] for f in fams}))
    df = pd.DataFrame(rows)
    fam_cols = [c for c in df.columns if c in set(fams) | set(K_UPTAKE) | set(SOLUTE)]
    df[fam_cols] = df[fam_cols].fillna(0.0)
    small = df['genome_bp'] < MIN_GENOME_BP
    print(f'excluded, assembly < {MIN_GENOME_BP:,} bp: {int(small.sum())}')
    df = df[~small].copy()
    df['n_k'] = (df[[c for c in K_UPTAKE if c in df]] > 0).sum(1)
    df['n_solute'] = (df[[c for c in SOLUTE if c in df]] > 0).sum(1)
    df['genus'] = df['genus'].fillna('unknown')
    return df, [c for c in fams if c in df]


def boot_diff(a: pd.DataFrame, b: pd.DataFrame, col: str, rng, n=4000):
    """Genus-cluster bootstrap CI for mean(a) - mean(b)."""
    def groups(d):
        return [g[col].to_numpy() for _, g in d.groupby('genus')]
    ga, gb = groups(a), groups(b)
    out = np.empty(n)
    for i in range(n):
        xa = np.concatenate([ga[j] for j in rng.integers(0, len(ga), len(ga))])
        xb = np.concatenate([gb[j] for j in rng.integers(0, len(gb), len(gb))])
        out[i] = xa.mean() - xb.mean()
    return np.percentile(out, [2.5, 97.5])


def compare(label, a, b, rng, out_rows):
    print(f'\n--- {label} ---')
    print(f'  salt-in n={len(a)} ({a["genus"].nunique()} genera)   '
          f'reference n={len(b)} ({b["genus"].nunique()} genera)')
    if len(a) < 3 or len(b) < 3:
        print('  too few strains')
        return
    print(f'  {"measure":<28}{"salt-in":>9}{"ref":>9}{"diff":>8}  95% CI (genus boot)   MWU p (genus means)')
    for col, name in (('n_k', 'K+ uptake (of 5)'), ('n_solute', 'compatible solutes (of 19)'),
                      ('n_pool', 'random pool (of 174)')):
        lo, hi = boot_diff(a, b, col, rng)
        ma, mb = a.groupby('genus')[col].mean(), b.groupby('genus')[col].mean()
        p = stats.mannwhitneyu(ma, mb, alternative='two-sided').pvalue
        print(f'  {name:<28}{a[col].mean():>9.2f}{b[col].mean():>9.2f}'
              f'{a[col].mean() - b[col].mean():>+8.2f}  [{lo:+.2f}, {hi:+.2f}]        {p:.2g}')
        out_rows.append(dict(comparison=label, measure=name, n_salt_in=len(a), n_ref=len(b),
                             genera_salt_in=a['genus'].nunique(), genera_ref=b['genus'].nunique(),
                             mean_salt_in=round(a[col].mean(), 3), mean_ref=round(b[col].mean(), 3),
                             diff=round(a[col].mean() - b[col].mean(), 3),
                             ci_lo=round(lo, 3), ci_hi=round(hi, 3), mwu_p_genus_means=p))
    d = pd.concat([a.assign(salt_in=1), b.assign(salt_in=0)])
    d['genome_mb'] = d['genome_bp'] / 1e6
    print('  adjusted for detectability: count ~ salt_in + random-pool count + genome Mb '
          '(OLS, genus-clustered SE)')
    for col, name in (('n_k', 'K+ uptake'), ('n_solute', 'compatible solutes')):
        try:
            m = smf.ols(f'{col} ~ salt_in + n_pool + genome_mb', d).fit(
                cov_type='cluster', cov_kwds={'groups': pd.factorize(d['genus'])[0]})
            lo, hi = m.conf_int().loc['salt_in']
            print(f'    {name:<20} salt_in coef {m.params["salt_in"]:+.2f} [{lo:+.2f}, {hi:+.2f}]  '
                  f'p={m.pvalues["salt_in"]:.2g}')
            out_rows.append(dict(comparison=label, measure=name + ' | pool, size',
                                 n_salt_in=len(a), n_ref=len(b),
                                 diff=round(m.params['salt_in'], 3), ci_lo=round(lo, 3),
                                 ci_hi=round(hi, 3), ols_p=m.pvalues['salt_in']))
        except Exception as e:  # singular design in a tiny stratum
            print(f'    {name}: model failed ({e})')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', default=None,
                    help='the BacDive data directory; or give its parts with the options below')
    ap.add_argument('--osmo-dir', default=None)
    ap.add_argument('--matrix', default=None, help='default: <osmo-dir>/feature_matrix.tsv')
    ap.add_argument('--classes', default=None, help='default: <dir>/bacdive_classes.tsv')
    ap.add_argument('--strains', default=None, help='default: <dir>/random_panel/strains')
    ap.add_argument('--control-pilot', default=None, help='default: <dir>/strategy_control/pilot_strains.tsv')
    ap.add_argument('--control-strains', default=None,
                    help='default: <dir>/strategy_control/random_panel/strains')
    ap.add_argument('--control-annotate', default=None,
                    help='default: <osmo-dir>/strategy_control/annotate')
    ap.add_argument('--out-dir', default=None, help='default: <osmo-dir>/strategy_control')
    ap.add_argument('--halophile-pct', type=float, default=15.0)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    args.osmo_dir = args.osmo_dir or args.dir
    if args.dir is None and not (args.matrix and args.classes and args.strains and args.control_pilot
                                 and args.control_strains and args.control_annotate and args.out_dir):
        ap.error('give --dir, or every one of --matrix, --classes, --strains, --control-* and --out-dir')
    rng = np.random.default_rng(args.seed)

    df, fams = load(args)
    arch = df[df['domain'] == 'Archaea']
    bact = df[df['domain'] == 'Bacteria']
    a_in, a_out = arch[arch['strategy'] == 'salt_in'], arch[arch['strategy'] == 'salt_out']
    b_in, b_out = bact[bact['strategy'] == 'salt_in'], bact[bact['strategy'] == 'salt_out']

    print('\n' + '=' * 76)
    print('WITHIN-DOMAIN SALT-IN CONTROL (post hoc extension, added 2026-11-14)')
    print('=' * 76)
    print('prediction: salt-in genomes carry K+ uptake families and few compatible-solute families')
    out_rows: list[dict] = []
    compare('Archaea: Halobacteria vs other archaea', a_in, a_out, rng, out_rows)
    print('\n  salt-out archaea by class (n, mean K+, mean solutes, mean pool):')
    for cls, g in a_out.groupby('tax_class'):
        print(f'    {str(cls):<40} n={len(g):<4} K {g["n_k"].mean():.2f}  '
              f'solutes {g["n_solute"].mean():.2f}  pool {g["n_pool"].mean():.0f}')
    compare('Bacteria: salt-in vs all salt-out bacteria', b_in, b_out, rng, out_rows)
    compare(f'Bacteria: salt-in vs salt-out with NaCl max >= {args.halophile_pct:g}%',
            b_in, b_out[b_out['nacl'] >= args.halophile_pct], rng, out_rows)
    hal = b_in[b_in['tax_class'] == 'Clostridia']
    if b_out['tax_class'].notna().any():
        compare('Bacteria: Halanaerobiales vs salt-out Clostridia',
                hal, b_out[b_out['tax_class'] == 'Clostridia'], rng, out_rows)
    print('\n  salt-in bacteria, per strain:')
    for r in b_in.sort_values('genus').itertuples():
        sol = [f for f in SOLUTE if getattr(r, f, 0) > 0]
        ks = [f for f in K_UPTAKE if getattr(r, f, 0) > 0]
        print(f'    {r.bacdive_id:>7} {r.genus:<18} K+: {",".join(ks) or "-":<24} solutes: {",".join(sol) or "-"}')

    # per-family prevalence by group
    groups = {'archaea_salt_in': a_in, 'archaea_salt_out': a_out,
              'bacteria_salt_in': b_in, 'bacteria_salt_out': b_out}
    prev = pd.DataFrame({name: (g[fams] > 0).mean() for name, g in groups.items()}).round(3)
    prev.insert(0, 'role', ['K_uptake' if f in K_UPTAKE else 'solute' if f in SOLUTE else 'other'
                            for f in prev.index])
    prev.index.name = 'family'
    print('\nper-family prevalence, K+ uptake and compatible-solute families:')
    print(prev[prev['role'] != 'other'].to_string())

    out_dir = args.out_dir or os.path.join(args.osmo_dir, 'strategy_control')
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(out_rows).to_csv(os.path.join(out_dir, 'strategy_groups.tsv'), sep='\t', index=False)
    prev.to_csv(os.path.join(out_dir, 'strategy_family_prevalence.tsv'), sep='\t')
    print(f'\nwrote {out_dir}/strategy_groups.tsv and strategy_family_prevalence.tsv')


if __name__ == '__main__':
    sys.exit(main())
