#!/usr/bin/env python3
"""Train the NaCl-max model on BacDive salt-out strains (Set A).

Steps, all written to --out:

  1. Training table     features.build_training_table() on --pv-dir;
                        training_table.tsv.
  2. Hyperparameters    random search over XGBoost settings, scored by
                        genus-blocked 5-fold CV RMSE on all of Set A;
                        hparam_search.tsv. The winner is used for steps 4-5
                        and for the final model.
  3. Nested CV          outer 5-fold genus-blocked (and family-blocked) CV
                        with the random search rerun inside each outer
                        training fold (4-fold genus-blocked). This is the
                        honest estimate for genomes from genera (families)
                        the model has never seen; oof_predictions.tsv.
  4. Ablations and      feature groups, target transform, genus weights,
     sensitivities      censored labels, plus ridge and the
                        HistGradientBoosting baseline of
                        analyze_extended.py. Fixed hyperparameters from
                        step 2 (so these are comparisons, not honest
                        estimates); cv_results.tsv.
  5. Final model        fitted on all of Set A; nacl_max.xgb.json.
  6. Uncertainty        split-conformal intervals: 10th/90th percentiles of
                        the nested genus-blocked OOF residuals, per bin of
                        the predicted value.
  7. Salt-in detector   logistic regression, salt-in vs salt-out on the same
                        features (all strains with a strategy, labelled or
                        not), genus-blocked OOF for its precision/recall.
                        Stored as plain coefficients in the schema.
  8. Schema             feature_schema.json: feature order, model settings,
                        CV metrics, intervals, applicability ranges, salt-in
                        detector, provenance.

Usage:
  python bin/nacl_model/train.py \\
      --pv-dir releases/v10/phenotype_validation \\
      --release releases/v10 --out releases/v10/nacl_model_build
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
from scipy import stats
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, RidgeCV
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import features as F  # noqa: E402

BINS = [(-np.inf, 1.0, 'non_tolerant'), (1.0, 3.0, 'slight'),
        (3.0, 15.0, 'moderate'), (15.0, np.inf, 'extreme')]
# predicted-value bins for the conformal intervals
INTERVAL_EDGES = [3.0, 6.0, 10.0]
ALPHA_LO, ALPHA_HI = 0.10, 0.90
# the 'baseline' model: genome size and GC, no osmo genes
BASELINE_FEATURES = ['total_proteins', 'gc']


def nacl_bin(v):
    for lo, hi, name in BINS:
        if lo < v <= hi:
            return name
    return BINS[0][2]


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def metrics(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    ss = ((y - y.mean()) ** 2).sum()
    return {
        'n': len(y),
        'spearman': stats.spearmanr(y, p)[0],
        'pearson': stats.pearsonr(y, p)[0],
        'mae': np.abs(y - p).mean(),
        'rmse': np.sqrt(((y - p) ** 2).mean()),
        'r2': 1 - ((y - p) ** 2).sum() / ss,
    }


# ---------------------------------------------------------------------------
# model wrappers
# ---------------------------------------------------------------------------

def sample_params(rng):
    return {
        'max_depth': int(rng.integers(2, 9)),
        'learning_rate': float(np.exp(rng.uniform(np.log(0.01), np.log(0.2)))),
        'n_estimators': int(rng.choice([150, 300, 500, 800, 1200])),
        'min_child_weight': float(np.exp(rng.uniform(np.log(1), np.log(50)))),
        'subsample': float(rng.uniform(0.5, 1.0)),
        'colsample_bytree': float(rng.uniform(0.3, 1.0)),
        'reg_lambda': float(np.exp(rng.uniform(np.log(0.1), np.log(30)))),
        'reg_alpha': float(rng.choice([0.0, 0.1, 1.0])),
    }


def fit_xgb(X, y, params, w=None, seed=0, threads=8):
    m = xgb.XGBRegressor(objective='reg:squarederror', tree_method='hist',
                         random_state=seed, n_jobs=threads, **params)
    m.fit(X, y, sample_weight=w)
    return m


def to_target(y, kind):
    return np.sqrt(y) if kind == 'sqrt' else y


def from_target(t, kind):
    t = np.asarray(t, float)
    if kind == 'sqrt':
        return np.clip(t, 0, None) ** 2
    return np.clip(t, 0, None)


def genus_weights(groups):
    n = pd.Series(groups).map(pd.Series(groups).value_counts()).to_numpy(float)
    w = 1 / np.sqrt(n)
    return w / w.mean()


def folds(groups, k, seed):
    return list(GroupKFold(n_splits=k, shuffle=True, random_state=seed).split(
        np.zeros(len(groups)), groups=groups))


def cv_xgb(X, y, groups, params, k=5, seed=0, target='raw', weights='none',
           threads=8, keep=None):
    """OOF predictions on the % scale. keep: optional bool mask of rows to
    train on (rows outside it are still predicted)."""
    oof = np.full(len(y), np.nan)
    for tr, te in folds(groups, k, seed):
        if keep is not None:
            tr = tr[keep[tr]]
        w = genus_weights(groups[tr]) if weights == 'genus' else None
        m = fit_xgb(X[tr], to_target(y[tr], target), params, w, seed, threads)
        oof[te] = from_target(m.predict(X[te]), target)
    return oof


def search(X, y, groups, n, k, seed, threads, label=''):
    rng = np.random.default_rng(seed)
    rows = []
    fl = folds(groups, k, seed)
    for i in range(n):
        p = sample_params(rng)
        oof = np.full(len(y), np.nan)
        for tr, te in fl:
            oof[te] = np.clip(fit_xgb(X[tr], y[tr], p, None, seed, threads).predict(X[te]), 0, None)
        mt = metrics(y, oof)
        rows.append({**p, **{f'cv_{a}': b for a, b in mt.items() if a != 'n'}})
        if label:
            log(f'  {label} config {i+1}/{n}: rmse {mt["rmse"]:.3f} rho {mt["spearman"]:.3f}')
    res = pd.DataFrame(rows).sort_values('cv_rmse').reset_index(drop=True)
    best = {c: res.loc[0, c] for c in sample_params(np.random.default_rng(0))}
    best = {c: (int(v) if c in ('max_depth', 'n_estimators') else float(v)) for c, v in best.items()}
    return best, res


# ---------------------------------------------------------------------------
# salt-in detector (plain logistic regression, shipped as coefficients)
# ---------------------------------------------------------------------------

def logit_features(X):
    return np.log1p(np.clip(X, 0, None))


def salt_in_detector(df, names, seed):
    d = df[df['strategy'].isin(['salt_in', 'salt_out'])].reset_index(drop=True)
    X = logit_features(d[names].to_numpy(float))
    y = (d['strategy'] == 'salt_in').to_numpy(int)
    g = d['genus'].to_numpy()

    def model():
        return make_pipeline(StandardScaler(), LogisticRegression(
            C=0.1, class_weight='balanced', max_iter=5000))

    oof = np.zeros(len(y))
    for tr, te in folds(g, 5, seed):
        oof[te] = model().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    rows = []
    for thr in [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]:
        pred = oof >= thr
        tp, fp = int((pred & (y == 1)).sum()), int((pred & (y == 0)).sum())
        rows.append({'threshold': thr, 'tp': tp, 'fp': fp, 'fn': int(((~pred) & (y == 1)).sum()),
                     'precision': tp / max(tp + fp, 1), 'recall': tp / max(y.sum(), 1)})
    perf = pd.DataFrame(rows)
    # threshold: the lowest with precision >= 0.8, else the highest-precision one
    ok = perf[perf['precision'] >= 0.8]
    thr = float(ok['threshold'].min()) if len(ok) else float(perf.sort_values('precision').iloc[-1]['threshold'])
    full = model().fit(X, y)
    sc, lr = full.named_steps['standardscaler'], full.named_steps['logisticregression']
    det = {
        'description': 'logistic regression on log1p(features), standardised; '
                       'p(salt-in) = sigmoid(intercept + sum(coef * (log1p(x) - mean) / scale))',
        'features': names,
        'mean': sc.mean_.tolist(), 'scale': sc.scale_.tolist(),
        'coef': lr.coef_[0].tolist(), 'intercept': float(lr.intercept_[0]),
        'threshold': thr,
        'n_salt_in': int(y.sum()), 'n_salt_out': int((y == 0).sum()),
        'cv_genus_blocked': perf.to_dict(orient='records'),
    }
    oof_df = d[['bacdive_id', 'species', 'genus', 'domain', 'strategy']].assign(p_salt_in=oof)
    return det, oof_df


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pv-dir', default=None, help='phenotype_validation dir of the release')
    ap.add_argument('--matrix', default=None, help='default: <pv-dir>/feature_matrix.tsv')
    ap.add_argument('--annotate', default=None, help='default: <pv-dir>/annotate')
    ap.add_argument('--strains', default=None, help='default: <pv-dir>/random_panel/strains')
    ap.add_argument('--release', required=True, help='release dir the annotate output was made with')
    ap.add_argument('--release-name', default=None, help="default: the --release directory's name")
    ap.add_argument('--osmotool-version', required=True,
                    help='osmotool version the training genomes were annotated with (recorded in the schema)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--n-search', type=int, default=40, help='random-search configs (step 2)')
    ap.add_argument('--n-search-inner', type=int, default=20, help='configs per outer fold (step 3)')
    ap.add_argument('--target', choices=['raw', 'sqrt'], default='raw')
    ap.add_argument('--weights', choices=['none', 'genus'], default='none')
    # Defaults from the ablations (README.md): systems and the breadth summary
    # add nothing over the family counts, and genome_mb nothing over
    # total_proteins -- and genome_bp is defined differently in the training
    # data (A/C/G/T only) and in osmotool's header (all characters).
    ap.add_argument('--drop', nargs='*', default=['genome_mb'],
                    help='features to leave out of the final model')
    ap.add_argument('--groups', nargs='*', default=['base', 'families'])
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--skip-nested', action='store_true', help='reuse oof_predictions.tsv')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rel_name = args.release_name or os.path.basename(os.path.normpath(args.release))

    # 1. table ---------------------------------------------------------------
    tpath = os.path.join(args.out, 'training_table.tsv')
    if os.path.exists(tpath):
        df = pd.read_csv(tpath, sep='\t', dtype={'bacdive_id': str})
        meta = json.load(open(os.path.join(args.out, 'training_table.meta.json')))
        families, systems = meta['families'], meta['systems']
        log(f'training table reused: {tpath}')
    else:
        log('building training table')
        if args.pv_dir is None and not (args.matrix and args.annotate and args.strains):
            sys.exit('give --pv-dir, or --matrix, --annotate and --strains')
        df, families, systems = F.build_training_table(args.pv_dir, matrix=args.matrix,
                                                       annotate=args.annotate, strains=args.strains)
        df.to_csv(tpath, sep='\t', index=False)
        json.dump({'families': families, 'systems': systems},
                  open(os.path.join(args.out, 'training_table.meta.json'), 'w'), indent=1)
    A = df[(df['strategy'] == 'salt_out') & df['nacl_max_pct'].notna()].reset_index(drop=True)
    # families never seen in Set A carry no information
    dead = [f for f in families if (A[f] > 0).sum() == 0]
    families = [f for f in families if f not in dead]
    all_names = F.feature_names(families, systems)
    names = [n for n in F.feature_names(families, systems, args.groups) if n not in args.drop]
    log(f'Set A: {len(A)} strains, {A["genus"].nunique()} genera, {A["family"].nunique()} families; '
        f'{len(names)} features (dropped all-zero: {dead}; --drop {args.drop})')

    y = A['nacl_max_pct'].to_numpy(float)
    gg = A['genus'].to_numpy()
    gf = A['family'].fillna('NA_' + A['genus']).to_numpy()
    cens = A['censored'].fillna(0).astype(int).to_numpy() == 1

    def Xof(cols):
        return A[cols].to_numpy(float)

    X = Xof(names)

    # 2. hyperparameter search -------------------------------------------------
    log(f'step 2: random search, {args.n_search} configs, genus-blocked 5-fold')
    best, sres = search(X, y, gg, args.n_search, 5, args.seed, args.threads, label='search')
    sres.to_csv(os.path.join(args.out, 'hparam_search.tsv'), sep='\t', index=False)
    log(f'best: {best}')

    # 3. nested CV ---------------------------------------------------------------
    oof_path = os.path.join(args.out, 'oof_predictions.tsv')
    def nested(Xm, label):
        """Nested CV on feature matrix Xm -> {block: (oof, fold_id)}"""
        out = {}
        for block, groups in [('genus', gg), ('family', gf)]:
            log(f'step 3: nested CV, {label}, {block}-blocked')
            oof = np.full(len(y), np.nan)
            fold_id = np.full(len(y), -1)
            for i, (tr, te) in enumerate(folds(groups, 5, args.seed)):
                p_in, _ = search(Xm[tr], y[tr], gg[tr], args.n_search_inner, 4, args.seed + i, args.threads)
                w = genus_weights(gg[tr]) if args.weights == 'genus' else None
                m = fit_xgb(Xm[tr], to_target(y[tr], args.target), p_in, w, args.seed, args.threads)
                oof[te] = from_target(m.predict(Xm[te]), args.target)
                fold_id[te] = i
                log(f'  outer fold {i+1}/5: rho {metrics(y[te], oof[te])["spearman"]:.3f} params {p_in}')
            out[block] = (oof, fold_id)
        return out

    nested_rows = []
    if args.skip_nested and os.path.exists(oof_path):
        oof_df = pd.read_csv(oof_path, sep='\t', dtype={'bacdive_id': str})
        log('step 3: reused nested OOF')
    else:
        oof_df = A[['bacdive_id', 'species', 'genus', 'family', 'nacl_max_pct',
                    'censored', 'nogrowth_at_pct']].copy()
        for block, (oof, fold_id) in nested(X, 'final features').items():
            oof_df[f'oof_{block}'] = oof
            oof_df[f'fold_{block}'] = fold_id
    # Baseline model: genome size and GC only, same nested procedure (own
    # hyperparameter search), to show how much the osmo genes add on top.
    if 'oof_baseline_genus' not in oof_df.columns:
        for block, (oof, _) in nested(Xof(BASELINE_FEATURES), 'baseline ' + '+'.join(BASELINE_FEATURES)).items():
            oof_df[f'oof_baseline_{block}'] = oof
    oof_df.to_csv(oof_path, sep='\t', index=False)
    for block in ['genus', 'family']:
        nested_rows.append({'experiment': 'nested_final', 'model': 'xgb', 'block': block,
                            'features': '+'.join(args.groups), 'n_features': len(names),
                            **metrics(y, oof_df[f'oof_{block}'])})
        log(f'nested {block}: {nested_rows[-1]}')
        nested_rows.append({'experiment': 'nested_baseline', 'model': 'xgb', 'block': block,
                            'features': '+'.join(BASELINE_FEATURES), 'n_features': len(BASELINE_FEATURES),
                            **metrics(y, oof_df[f'oof_baseline_{block}'])})
        log(f'nested baseline {block}: {nested_rows[-1]}')

    # 4. ablations / sensitivities / baselines (fixed params) --------------------
    log('step 4: ablations')
    exps = {
        'base': F.feature_names(families, systems, ['base']),
        'base+families': F.feature_names(families, systems, ['base', 'families']),
        'base+families+summary': F.feature_names(families, systems, ['base', 'families', 'summary']),
        'full': all_names,
        'full-genome_mb': [n for n in all_names if n != 'genome_mb'],
        'full-summary': [n for n in all_names if n not in ('n_families_present', 'total_osmo_copies')],
        'full-systems': F.feature_names(families, systems, ['base', 'families', 'summary']),
        'families+summary+systems (no base)': F.feature_names(families, systems, ['families', 'summary', 'systems']),
        'selected': names,
    }
    rows = []
    for block, groups in [('genus', gg), ('family', gf)]:
        for ename, cols in exps.items():
            oof = cv_xgb(Xof(cols), y, groups, best, seed=args.seed, threads=args.threads)
            rows.append({'experiment': 'ablation', 'model': 'xgb', 'block': block,
                         'features': ename, 'n_features': len(cols), **metrics(y, oof)})
        for tname in ['sqrt']:
            oof = cv_xgb(X, y, groups, best, seed=args.seed, target=tname, threads=args.threads)
            rows.append({'experiment': 'target=sqrt', 'model': 'xgb', 'block': block,
                         'features': 'selected', 'n_features': len(names), **metrics(y, oof)})
        oof = cv_xgb(X, y, groups, best, seed=args.seed, weights='genus', threads=args.threads)
        rows.append({'experiment': 'weights=1/sqrt(n_genus)', 'model': 'xgb', 'block': block,
                     'features': 'selected', 'n_features': len(names), **metrics(y, oof)})
        oof = cv_xgb(X, y, groups, best, seed=args.seed, threads=args.threads, keep=~cens)
        rows.append({'experiment': 'train without censored', 'model': 'xgb', 'block': block,
                     'features': 'selected', 'n_features': len(names), **metrics(y, oof)})
        rows.append({'experiment': 'eval uncensored only', 'model': 'xgb', 'block': block,
                     'features': 'selected', 'n_features': len(names),
                     **metrics(y[~cens], oof_df[f'oof_{block}'].to_numpy()[~cens])})
        # baselines: ridge on log1p, HGB as in analyze_extended.py
        Xl = np.log1p(np.clip(X, 0, None))
        oof_r, oof_h = np.full(len(y), np.nan), np.full(len(y), np.nan)
        for tr, te in folds(groups, 5, args.seed):
            r = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25))).fit(Xl[tr], y[tr])
            oof_r[te] = np.clip(r.predict(Xl[te]), 0, None)
            h = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05,
                                              random_state=args.seed).fit(Xl[tr], y[tr])
            oof_h[te] = np.clip(h.predict(Xl[te]), 0, None)
        rows.append({'experiment': 'baseline', 'model': 'ridge(log1p)', 'block': block,
                     'features': 'selected', 'n_features': len(names), **metrics(y, oof_r)})
        rows.append({'experiment': 'baseline', 'model': 'hgb(log1p, default)', 'block': block,
                     'features': 'selected', 'n_features': len(names), **metrics(y, oof_h)})
        rows.append({'experiment': 'baseline', 'model': 'predict global mean', 'block': block,
                     'features': '-', 'n_features': 0, **metrics(y, np.full(len(y), y.mean()))})
        log(f'  {block} done')
    cv = pd.DataFrame(nested_rows + rows)
    cv.to_csv(os.path.join(args.out, 'cv_results.tsv'), sep='\t', index=False, float_format='%.4f')
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print(cv[['experiment', 'model', 'block', 'features', 'spearman', 'mae', 'rmse', 'r2']]
              .to_string(index=False, float_format='%.3f'))

    # per true-bin metrics of the nested genus OOF
    bin_rows = []
    for block in ['genus', 'family']:
        o = oof_df[f'oof_{block}'].to_numpy()
        tb = np.array([nacl_bin(v) for v in y])
        pb = np.array([nacl_bin(v) for v in o])
        for _, _, b in BINS:
            s = tb == b
            if s.sum():
                bin_rows.append({'block': block, 'true_bin': b, 'n': int(s.sum()),
                                 'mean_true': y[s].mean(), 'mean_pred': o[s].mean(),
                                 'mae': np.abs(y[s] - o[s]).mean(),
                                 'bin_accuracy': (pb[s] == b).mean()})
    pd.DataFrame(bin_rows).to_csv(os.path.join(args.out, 'cv_by_bin.tsv'), sep='\t',
                                  index=False, float_format='%.3f')

    # 5. final model -----------------------------------------------------------
    log('step 5: final model')
    w = genus_weights(gg) if args.weights == 'genus' else None
    final = fit_xgb(X, to_target(y, args.target), best, w, args.seed, args.threads)
    booster = final.get_booster()
    booster.feature_names = names
    mpath = os.path.join(args.out, 'nacl_max.xgb.json')
    booster.save_model(mpath)

    # 6. conformal intervals from nested genus OOF -------------------------------
    o = oof_df['oof_genus'].to_numpy()
    res = y - o
    edges = [-np.inf] + INTERVAL_EDGES + [np.inf]
    intervals = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (o > lo) & (o <= hi)
        intervals.append({'pred_gt': None if np.isinf(lo) else lo,
                          'pred_le': None if np.isinf(hi) else hi, 'n': int(s.sum()),
                          'resid_lo': float(np.quantile(res[s], ALPHA_LO)),
                          'resid_hi': float(np.quantile(res[s], ALPHA_HI))})

    def apply_iv(pred):
        lo_, hi_ = np.zeros(len(pred)), np.zeros(len(pred))
        for iv in intervals:
            a = -np.inf if iv['pred_gt'] is None else iv['pred_gt']
            b = np.inf if iv['pred_le'] is None else iv['pred_le']
            s = (pred > a) & (pred <= b)
            lo_[s] = np.clip(pred[s] + iv['resid_lo'], 0, None)
            hi_[s] = pred[s] + iv['resid_hi']
        return lo_, hi_
    lo_g, hi_g = apply_iv(o)
    of = oof_df['oof_family'].to_numpy()
    lo_f, hi_f = apply_iv(of)
    coverage = {'genus_oof_in_sample': float(((y >= lo_g) & (y <= hi_g)).mean()),
                'family_oof_with_genus_intervals': float(((y >= lo_f) & (y <= hi_f)).mean()),
                'nominal': ALPHA_HI - ALPHA_LO}
    log(f'interval coverage: {coverage}')

    # 7. salt-in detector ---------------------------------------------------------
    log('step 7: salt-in detector')
    det_names = [n for n in all_names if n not in ('n_contigs', 'genome_mb')]
    det, det_oof = salt_in_detector(df, det_names, args.seed)
    det_oof.to_csv(os.path.join(args.out, 'salt_in_detector_oof.tsv'), sep='\t', index=False)
    log(f'salt-in threshold {det["threshold"]}: '
        f'{[r for r in det["cv_genus_blocked"] if r["threshold"] == det["threshold"]]}')

    # salt-in strains with a label: out-of-domain readout of the final model
    S = df[(df['strategy'] == 'salt_in') & df['nacl_max_pct'].notna()]
    if len(S):
        ps = from_target(final.predict(S[names].to_numpy(float)), args.target)
        S = S[['bacdive_id', 'species', 'genus', 'domain', 'nacl_max_pct']].assign(pred=ps)
        S.to_csv(os.path.join(args.out, 'salt_in_predictions.tsv'), sep='\t', index=False,
                 float_format='%.2f')
        log(f'salt-in (n={len(S)}): median true {S.nacl_max_pct.median():.1f}, '
            f'median pred {S.pred.median():.1f}')

    # 8. schema --------------------------------------------------------------------
    q = A[names].quantile([0.01, 0.99])
    sha = hashlib.sha256(open(tpath, 'rb').read()).hexdigest()
    mt_g = metrics(y, oof_df['oof_genus'])
    mt_f = metrics(y, oof_df['oof_family'])
    schema = {
        'model': 'nacl_max', 'schema_version': 1,
        'description': 'Predicts the highest NaCl concentration (% w/v) at which a salt-out '
                       'prokaryote grows, from osmotool annotate output. Trained on BacDive.',
        'model_file': 'nacl_max.xgb.json',
        'created': time.strftime('%Y-%m-%d'),
        'refdb_release': rel_name,
        'annotate': {'method': 'hmm', 'prodigal': '-p meta', 'osmotool_version_training': args.osmotool_version,
                     'min_osmotool_version_inference': '0.8.0'},
        'features': names,
        'families': families, 'systems': systems,
        'feature_definitions': {
            'total_proteins': "gene_counts.tsv header '# total_proteins'",
            'gc': "gene_counts.tsv header '# gc'",
            'n_contigs': "gene_counts.tsv header '# n_contigs'",
            'genome_mb': "gene_counts.tsv header '# genome_bp' / 1e6",
            '<family>': 'gene_counts.tsv raw_count; 0 when the family is absent',
            'n_families_present': 'number of schema families with raw_count > 0',
            'total_osmo_copies': 'sum of raw_count over schema families',
            'sys_<s>_present': 'systems.tsv present (True=1)',
            'sys_<s>_subunits': 'systems.tsv subunit_count',
        },
        'target': {'name': 'nacl_max_pct', 'unit': '% (w/v) NaCl', 'transform': args.target,
                   'sample_weights': args.weights,
                   'clip_min': 0.0},
        'xgboost_params': {'objective': 'reg:squarederror', 'tree_method': 'hist', **best},
        'training': {'n_strains': int(len(A)), 'n_genera': int(A['genus'].nunique()),
                     'n_families': int(A['family'].nunique()),
                     'set': "BacDive strains with a quantitative NaCl max, strategy salt_out",
                     'training_table_sha256': sha, 'target_mean': float(y.mean()),
                     'target_quantiles': {str(k): float(v) for k, v in
                                          zip([0.1, 0.25, 0.5, 0.75, 0.9],
                                              np.quantile(y, [0.1, 0.25, 0.5, 0.75, 0.9]))}},
        'cv': {'genus_blocked': mt_g, 'family_blocked': mt_f,
               'baseline_features': BASELINE_FEATURES,
               'baseline_genus_blocked': metrics(y, oof_df['oof_baseline_genus']),
               'baseline_family_blocked': metrics(y, oof_df['oof_baseline_family']),
               'note': 'nested 5-fold CV, hyperparameters re-searched inside each outer fold'},
        'bins': [{'name': b, 'gt': None if np.isinf(lo) else lo, 'le': None if np.isinf(hi) else hi}
                 for lo, hi, b in BINS],
        'intervals': {'method': 'split-conformal, residual quantiles of nested genus-blocked OOF '
                                'predictions per predicted-value bin',
                      'alpha_lo': ALPHA_LO, 'alpha_hi': ALPHA_HI, 'bins': intervals,
                      'coverage': coverage},
        'applicability': {
            'feature_q01': {n: float(q.loc[0.01, n]) for n in names},
            'feature_q99': {n: float(q.loc[0.99, n]) for n in names},
            'max_features_outside': 3,
            'min_total_proteins': float(A['total_proteins'].min()),
            'n_contigs_warn': float(A['n_contigs'].quantile(0.99)),
        },
        'salt_in_detector': det,
        'software': {'python': platform.python_version(), 'xgboost': xgb.__version__,
                     'scikit-learn': sklearn.__version__, 'numpy': np.__version__,
                     'pandas': pd.__version__},
    }
    json.dump(schema, open(os.path.join(args.out, 'feature_schema.json'), 'w'), indent=1,
              default=float)
    log(f'wrote {mpath} and feature_schema.json')


if __name__ == '__main__':
    main()
