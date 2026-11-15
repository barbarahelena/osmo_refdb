#!/usr/bin/env python3
"""The parts of the prespecified plan that analyze.py does not cover.

analyze.py is the dependency-free readout written for the pilot. This adds
what the plan in README.md commits to but analyze.py skips:

  1. Baselines   genome size and GC alongside total_proteins, and the panel's
                 correlation with all three partialled out jointly.
  2. Random      where the osmo panel falls among many random panels of
     panels      non-osmotic Pfam families (pool from build_random_pool.py),
                 drawn unconstrained and prevalence-matched. Same size as
                 the panel's breadth measure: 44 columns.
  3. Set B       nacl_max ~ panel + covariates + (1|genus) mixed model, and the
     model       within-genus (genus-demeaned) estimate, both also run on the
                 random panels so the real panel has a null to sit in.
  4. Secondary   all strains incl. salt-in, with a strategy indicator and
                 panel x strategy interaction.
  5. Grouped CV  out-of-fold prediction with folds blocked by genus and by
                 taxonomic family: baselines alone vs baselines + panel vs
                 baselines + random panels.
  6. Family      per-family association, raw, baseline-partialled and
     ranking     within-genus.
  7. Salt-in     whether the solute "absence" in salt-in strains is specific,
     check       or HMM detection being generally weaker in those genomes.

Needs numpy, pandas, scikit-learn and statsmodels (not in the image).

Usage:
  python analyze_extended.py --dir releases/v8/phenotype_validation_full
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import statsmodels.formula.api as smf

K_UPTAKE = ['trkA', 'trkH', 'ktrA', 'ktrB', 'ktrD']
SOLUTE = ['ectA', 'ectB', 'ectC', 'betL', 'betA', 'betB', 'otsA', 'otsB',
          'proX', 'proP', 'opuAA', 'opuAB', 'opuAC', 'opuBA', 'opuBB', 'opuBC',
          'opuCA', 'opuCB', 'opuCC']
META = {'bacdive_id', 'species', 'genus', 'family', 'domain', 'strategy',
        'nacl_max_pct', 'censored', 'nogrowth_at_pct', 'asm_level',
        'total_proteins', 'n_families_present', 'total_osmo_copies'}
BASE = ['log_proteins', 'genome_mb', 'gc']


def hdr(title):
    print()
    print('=' * 76)
    print(title)
    print('=' * 76)


def rank(a):
    return stats.rankdata(a, axis=0)


def resid(y, X):
    """Residuals of y (n,) or (n,k) after OLS on X (n,p) plus intercept."""
    X1 = np.column_stack([np.ones(len(X)), X])
    beta, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return y - X1 @ beta


def corr_cols(A, b):
    """Pearson r between each column of A (n,k) and vector b (n,)."""
    A = A - A.mean(0)
    b = b - b.mean()
    den = np.sqrt((A * A).sum(0) * (b * b).sum())
    with np.errstate(invalid='ignore', divide='ignore'):
        return (A * b[:, None]).sum(0) / den


def partial_spearman(X, y, Z):
    """Rank-level partial correlation of each column of X with y, given Z."""
    rX, ry, rZ = rank(X), rank(y), rank(Z)
    return corr_cols(resid(rX, rZ), resid(ry, rZ))


def demean(v, groups):
    """Subtract group means (v may be (n,) or (n,k))."""
    df = pd.DataFrame(v)
    return (df - df.groupby(groups.values).transform('mean')).to_numpy().squeeze()


def within_partial(X, y, Z, groups):
    """Within-genus partial correlation: genus-demean everything, then partial."""
    Xd, yd, Zd = demean(X, groups), demean(y, groups), demean(Z, groups)
    if Xd.ndim == 1:
        Xd = Xd[:, None]
    return corr_cols(resid(Xd, Zd), resid(yd, Zd))


def within_genus_median_rho(B, y, groups):
    """Median over genera of Spearman(breadth, y); B is (n,k) for k panels."""
    df = pd.DataFrame(B)
    g = groups.values
    rB = df.groupby(g).rank().to_numpy()
    ry = pd.Series(y).groupby(g).rank().to_numpy()
    rB = rB - pd.DataFrame(rB).groupby(g).transform('mean').to_numpy()
    ry = ry - pd.Series(ry).groupby(g).transform('mean').to_numpy()
    num = pd.DataFrame(rB * ry[:, None]).groupby(g).sum()
    sb = pd.DataFrame(rB * rB).groupby(g).sum()
    sy = pd.Series(ry * ry).groupby(g).sum()
    with np.errstate(invalid='ignore', divide='ignore'):
        r = num.to_numpy() / np.sqrt(sb.to_numpy() * sy.to_numpy()[:, None])
    return np.nanmedian(r, axis=0)


def cluster_boot(fn, n, groups, reps, rng):
    """95% CI of fn(idx) resampling whole groups."""
    by = pd.Series(np.arange(n)).groupby(groups.values).apply(np.array).tolist()
    out = []
    for _ in range(reps):
        pick = rng.integers(0, len(by), len(by))
        out.append(fn(np.concatenate([by[i] for i in pick])))
    return np.nanpercentile(out, [2.5, 97.5])


def pct_rank(real, null):
    """Share of null panels scoring >= the real panel (one-sided, +1 smoothed)."""
    null = np.asarray(null)
    null = null[~np.isnan(null)]
    return (1 + (null >= real).sum()) / (1 + len(null))


def load(args):
    fm = pd.read_csv(os.path.join(args.dir, 'feature_matrix.tsv'), sep='\t',
                     dtype={'bacdive_id': str})
    osmo = [c for c in fm.columns if c not in META]

    rp = os.path.join(args.dir, 'random_panel')
    pool = pd.read_csv(os.path.join(rp, 'pool_families.tsv'), sep='\t')
    stats_rows, count_rows = [], {}
    for d in glob.glob(os.path.join(rp, 'strains', '*')):
        bid = os.path.basename(d)
        s, c = os.path.join(d, f'{bid}.genome_stats.tsv'), os.path.join(d, f'{bid}.pool_counts.tsv')
        if not (os.path.exists(s) and os.path.exists(c)):
            continue
        r = pd.read_csv(s, sep='\t').iloc[0].to_dict()
        r['bacdive_id'] = bid
        stats_rows.append(r)
        cc = pd.read_csv(c, sep='\t')
        count_rows[bid] = dict(zip(cc['family'], cc['n_proteins']))
    gs = pd.DataFrame(stats_rows).rename(columns={'total_proteins': 'prodigal_proteins'})
    rc = (pd.DataFrame.from_dict(count_rows, orient='index')
          .reindex(columns=pool['name']).fillna(0))
    rc.index.name = 'bacdive_id'
    rc.columns = [f'R_{c}' for c in rc.columns]

    df = fm.merge(gs, on='bacdive_id', how='inner').merge(
        rc.reset_index(), on='bacdive_id', how='inner')
    df = df[df['nacl_max_pct'].notna()].copy()
    df['log_proteins'] = np.log(df['total_proteins'])
    df['genome_mb'] = df['genome_bp'] / 1e6
    df['salt_in'] = (df['strategy'] == 'salt_in').astype(int)
    rand = list(rc.columns)
    mismatch = (df['prodigal_proteins'] != df['total_proteins']).sum()
    return df, osmo, rand, len(fm), mismatch


def set_b(df, args):
    g = df.groupby('genus')['nacl_max_pct']
    ok = g.transform('size').ge(args.min_genus_n) & \
        (g.transform('max') - g.transform('min')).ge(args.min_genus_spread)
    return df[ok].copy()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--n-panels', type=int, default=1000)
    ap.add_argument('--boot', type=int, default=2000)
    ap.add_argument('--cv-draws', type=int, default=20)
    ap.add_argument('--min-genus-n', type=int, default=3)
    ap.add_argument('--min-genus-spread', type=float, default=4.0)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    df_all, osmo, rand, n_fm, mismatch = load(args)
    A = df_all[df_all['salt_in'] == 0].reset_index(drop=True)
    B = set_b(A, args).reset_index(drop=True)
    y = A['nacl_max_pct'].to_numpy(float)
    Z = A[BASE].to_numpy(float)
    gen = A['genus']
    P = len(osmo)

    hdr('DATA')
    print(f'feature matrix strains          : {n_fm}')
    print(f'with random-panel + genome stats: {len(df_all)} (quantitative NaCl max)')
    print(f'Set A (salt-out)                : {len(A)} strains, {gen.nunique()} genera')
    print(f'Set B (within-genus contrast)   : {len(B)} strains, {B["genus"].nunique()} genera')
    print(f'salt-in                         : {int(df_all["salt_in"].sum())}')
    print(f'osmo panel columns              : {P}    random pool families: {len(rand)}')
    print(f'Prodigal protein count != osmotool total_proteins: {mismatch} strains '
          '(should be 0: same Prodigal, same mode)')

    # ---------------------------------------------------------------- 1
    hdr('1. BASELINES -- Set A, Spearman rho vs NaCl max [95% genus-cluster boot]')
    breadth = A['n_families_present'].to_numpy(float)

    def show(label, fn):
        est = fn(np.arange(len(A)))
        lo, hi = cluster_boot(fn, len(A), gen, args.boot, rng)
        print(f'  {label:<44} {est:+.3f}  [{lo:+.3f}, {hi:+.3f}]')

    for col, label in (('total_proteins', 'total_proteins'), ('genome_mb', 'genome size'),
                       ('gc', 'GC'), ('n_contigs', 'n_contigs (assembly quality; diagnostic)')):
        v = A[col].to_numpy(float)
        show(label, lambda i, v=v: stats.spearmanr(v[i], y[i])[0])
    show('panel breadth', lambda i: stats.spearmanr(breadth[i], y[i])[0])
    show('panel breadth | proteins + genome size + GC',
         lambda i: partial_spearman(breadth[i][:, None], y[i], Z[i])[0])

    # ---------------------------------------------------------------- 2
    hdr(f'2. RANDOM PANELS -- {args.n_panels} panels of {P} non-osmotic Pfam families')
    Rpres = (A[rand].to_numpy() > 0)
    prev_r = Rpres.mean(0)
    prev_o = (A[osmo].to_numpy() > 0).mean(0)
    print(f'  prevalence, osmo panel   : median {np.median(prev_o):.2f} '
          f'(IQR {np.percentile(prev_o, 25):.2f}-{np.percentile(prev_o, 75):.2f})')
    print(f'  prevalence, random pool  : median {np.median(prev_r):.2f} '
          f'(IQR {np.percentile(prev_r, 25):.2f}-{np.percentile(prev_r, 75):.2f})')

    def draw_unconstrained():
        return rng.choice(len(rand), P, replace=False)

    order = np.argsort(prev_r)

    def draw_matched(k=5):
        """Per osmo family, a random pick among the k unused pool families
        nearest in prevalence."""
        used, pick = set(), []
        for p in rng.permutation(prev_o):
            cand = [j for j in order[np.argsort(np.abs(prev_r[order] - p))] if j not in used][:k]
            j = rng.choice(cand)
            used.add(j)
            pick.append(j)
        return np.array(pick)

    yB = B['nacl_max_pct'].to_numpy(float)
    ZB = B[BASE].to_numpy(float)
    genB = B['genus']
    RpresB = (B[rand].to_numpy() > 0)
    breadthB = B['n_families_present'].to_numpy(float)

    real = dict(
        rho=stats.spearmanr(breadth, y)[0],
        partial=partial_spearman(breadth[:, None], y, Z)[0],
        within=within_partial(breadthB, yB, ZB, genB)[0],
        median_within=within_genus_median_rho(breadthB[:, None], yB, genB)[0],
    )
    null_rows = []
    for kind, draw in (('unconstrained', draw_unconstrained), ('prevalence-matched', draw_matched)):
        picks = [draw() for _ in range(args.n_panels)]
        BA = np.column_stack([Rpres[:, p].sum(1) for p in picks]).astype(float)
        BB = np.column_stack([RpresB[:, p].sum(1) for p in picks]).astype(float)
        res = dict(
            rho=corr_cols(rank(BA), rank(y)),
            partial=partial_spearman(BA, y, Z),
            within=within_partial(BB, yB, ZB, genB),
            median_within=within_genus_median_rho(BB, yB, genB),
        )
        print(f'\n  {kind}')
        print(f'  {"statistic":<38}{"osmo":>8}{"null median":>13}{"null 95th":>11}{"p(null>=osmo)":>15}')
        for key, label in (('rho', 'Set A rho'), ('partial', 'Set A rho | baselines'),
                           ('within', 'Set B within-genus partial r'),
                           ('median_within', 'Set B median within-genus rho')):
            nl = res[key]
            print(f'  {label:<38}{real[key]:>+8.3f}{np.nanmedian(nl):>+13.3f}'
                  f'{np.nanpercentile(nl, 95):>+11.3f}{pct_rank(real[key], nl):>15.4f}')
        for i in range(args.n_panels):
            null_rows.append({'kind': kind, 'panel': i, **{k: res[k][i] for k in res}})
    pd.DataFrame(null_rows).to_csv(os.path.join(args.dir, 'null_panels.tsv'), sep='\t', index=False)
    print('\n  p is one-sided: the share of random panels scoring at least as high as')
    print('  the osmo panel (+1 smoothed, so the floor is 1/(n_panels+1)).')
    print('  Pool panels overlap (44 drawn from a fixed pool), so the null is somewhat')
    print('  narrower than independent draws from all of Pfam would give.')

    # ---------------------------------------------------------------- 3
    hdr('3. SET B MIXED MODEL -- nacl_max ~ panel + baselines + (1|genus)')
    for name, d in (('Set B', B), ('Set A', A)):
        d = d.copy()
        for c in ['n_families_present'] + BASE:
            d[c + '_z'] = (d[c] - d[c].mean()) / d[c].std()
        m = smf.mixedlm('nacl_max_pct ~ n_families_present_z + log_proteins_z + genome_mb_z + gc_z',
                        d, groups=d['genus']).fit(reml=True, method='lbfgs')
        b, se = m.params['n_families_present_z'], m.bse['n_families_present_z']
        print(f'  {name}: panel breadth, per SD = {b:+.3f} % NaCl  '
              f'[{b-1.96*se:+.3f}, {b+1.96*se:+.3f}]  p={m.pvalues["n_families_present_z"]:.1e}'
              f'  (n={len(d)}, genera={d["genus"].nunique()}, '
              f'genus variance={float(m.cov_re.iloc[0, 0]):.2f}, residual={m.scale:.2f})')
    print('  A random-intercept model blends within- and between-genus variation;')
    print('  the genus-demeaned partial r in section 2 is the purely within-genus')
    print('  estimate, and is the one with a random-panel null.')

    # ---------------------------------------------------------------- 4
    hdr('4. SECONDARY -- all strains, strategy indicator + panel x strategy')
    d = df_all.copy()
    for c in ['n_families_present'] + BASE:
        d[c + '_z'] = (d[c] - d[c].mean()) / d[c].std()
    m = smf.mixedlm('nacl_max_pct ~ n_families_present_z * salt_in + log_proteins_z + genome_mb_z + gc_z',
                    d, groups=d['genus']).fit(reml=True, method='lbfgs')
    for term in ('n_families_present_z', 'salt_in', 'n_families_present_z:salt_in'):
        b, se = m.params[term], m.bse[term]
        print(f'  {term:<30} {b:+.3f}  [{b-1.96*se:+.3f}, {b+1.96*se:+.3f}]  p={m.pvalues[term]:.1e}')
    n_si_gen = d.loc[d['salt_in'] == 1, 'genus'].nunique()
    shared = set(d.loc[d['salt_in'] == 1, 'genus']) & set(d.loc[d['salt_in'] == 0, 'genus'])
    print(f'  salt-in: {int(d["salt_in"].sum())} strains in {n_si_gen} genera, '
          f'{len(shared)} genera shared with salt-out -- strategy is almost entirely')
    print('  between-genus, so the salt_in terms are carried by few clusters.')

    # ---------------------------------------------------------------- 5
    hdr('5. GROUPED CROSS-VALIDATION -- out-of-fold Spearman (R^2), 5 folds')

    def folds(groups, k=5):
        u = groups.unique()
        assign = dict(zip(rng.permutation(u), np.arange(len(u)) % k))
        return groups.map(assign).to_numpy()

    def oof(X, model):
        pred = np.empty(len(y))
        for f in range(5):
            tr, te = fold != f, fold == f
            pred[te] = model().fit(X[tr], y[tr]).predict(X[te])
        return stats.spearmanr(pred, y)[0], 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()

    ridge = lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25)))
    hgb = lambda: HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05,
                                                random_state=args.seed)
    Xo = np.log1p(A[osmo].to_numpy(float))
    Xr_all = np.log1p(A[rand].to_numpy(float))
    rand_draws = [draw_matched() for _ in range(args.cv_draws)]
    cv_rows = []
    for block in ('genus', 'family'):
        fold = folds(A[block])
        print(f'\n  blocked by {block} ({A[block].nunique()} groups)')
        for mname, model in (('ridge', ridge), ('gradient boosting', hgb)):
            base = oof(Z, model)
            panel = oof(np.column_stack([Z, Xo]), model)
            rnd = np.array([oof(np.column_stack([Z, Xr_all[:, p]]), model) for p in rand_draws])
            print(f'    {mname}')
            print(f'      baselines only                 {base[0]:+.3f} ({base[1]:+.3f})')
            print(f'      baselines + osmo panel         {panel[0]:+.3f} ({panel[1]:+.3f})')
            print(f'      baselines + random panel       {np.median(rnd[:, 0]):+.3f} ({np.median(rnd[:, 1]):+.3f})'
                  f'   median of {args.cv_draws} prevalence-matched draws; '
                  f'range {rnd[:, 0].min():+.3f} to {rnd[:, 0].max():+.3f}')
            cv_rows += [dict(block=block, model=mname, features='baselines', spearman=base[0], r2=base[1]),
                        dict(block=block, model=mname, features='baselines+osmo', spearman=panel[0], r2=panel[1])]
            cv_rows += [dict(block=block, model=mname, features=f'baselines+random_{i}',
                             spearman=s, r2=r) for i, (s, r) in enumerate(rnd)]
    pd.DataFrame(cv_rows).to_csv(os.path.join(args.dir, 'cv_results.tsv'), sep='\t', index=False)

    # ---------------------------------------------------------------- 6
    hdr('6. FAMILY RANKING -- presence vs NaCl max, Set A / Set B')
    pres = (A[osmo].to_numpy() > 0).astype(float)
    presB = (B[osmo].to_numpy() > 0).astype(float)
    raw = corr_cols(rank(pres), rank(y))
    part = partial_spearman(pres, y, Z)
    within = within_partial(presB, yB, ZB, genB)
    rk = pd.DataFrame({'family': osmo, 'prevalence': pres.mean(0), 'rho': raw,
                       'rho_given_baselines': part, 'within_genus_r': within})
    rk['group'] = np.where(rk['family'].isin(K_UPTAKE), 'K+ uptake',
                           np.where(rk['family'].isin(SOLUTE), 'solute', ''))
    rk = rk.sort_values('within_genus_r', ascending=False)
    rk.to_csv(os.path.join(args.dir, 'family_ranking.tsv'), sep='\t', index=False)
    print(f'  {"family":<18}{"prev":>6}{"rho":>8}{"|base":>8}{"within":>8}  group')
    for _, r in rk.iterrows():
        print(f'  {r.family:<18}{r.prevalence:>6.2f}{r.rho:>+8.3f}{r.rho_given_baselines:>+8.3f}'
              f'{r.within_genus_r:>+8.3f}  {r.group}')
    print('  Sorted by within-genus r. mrpA-G are subunits of one antiporter complex and')
    print('  co-occur by construction: several mrp rows are one observation, not several.')

    # ---------------------------------------------------------------- 7
    hdr('7. SALT-IN CHECK -- is solute absence specific, or weaker detection overall?')
    si = df_all[df_all['salt_in'] == 1]
    arch_out = df_all[(df_all['salt_in'] == 0) & (df_all['domain'] == 'Archaea')]

    def frac(d, cols):
        return (d[cols].to_numpy() > 0).mean() if len(d) else float('nan')
    print(f'  {"fraction present":<30}{"salt-in":>10}{"salt-out":>10}{"salt-out archaea":>18}')
    for label, cols in (('K+ uptake', K_UPTAKE), ('compatible solutes', SOLUTE),
                        ('random pool families', rand)):
        print(f'  {label:<30}{frac(si, cols):>10.2f}{frac(A, cols):>10.2f}{frac(arch_out, cols):>18.2f}')
    print(f'  (salt-in n={len(si)}, salt-out n={len(A)}, salt-out archaea n={len(arch_out)})')
    print('  If random families drop in salt-in strains as much as solutes do, the')
    print('  solute result is not specific to the solute machinery.')


if __name__ == '__main__':
    main()
