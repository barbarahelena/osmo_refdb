#!/usr/bin/env python3
"""What the NaCl-max model weighs, and how stable that is.

Reads the output of train.py (--build) and writes to <build>/explain/:

  shap_importance.tsv      per feature: mean |SHAP| (% NaCl), direction
                           (Spearman of feature value vs SHAP), prevalence,
                           XGBoost gain / cover / weight importance
  group_importance.tsv     same, with SHAP summed over functional groups
  importance_stability.tsv mean |SHAP| per genus-blocked CV fold (refit with
                           the final hyperparameters), rank per fold
  interactions.tsv         top feature pairs by mean |SHAP interaction|
  ridge_coefficients.tsv   standardised ridge on log1p features (a linear
                           reference), next to the within-genus correlations
                           of phenotype_validation/family_ranking.tsv
  examples.tsv             per-strain SHAP breakdown for the example strains
  *.png                    bar, beeswarm, dependence, group, interaction and
                           waterfall plots
  README.md                plain-language summary of all the above (shipped as
                           <release>/nacl_model/README.md)

SHAP values are XGBoost's own TreeSHAP (Booster.predict(pred_contribs=True)):
for every strain the bias (mean training prediction) plus the per-feature
contributions add up exactly to the prediction, in % NaCl.

Usage:
  python bin/nacl_model/explain.py --build releases/v10/nacl_model_build \\
      --pv-dir releases/v10/phenotype_validation
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xgboost as xgb  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.model_selection import GroupKFold  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

GROUPS = {
    'genome size / composition': ['total_proteins', 'gc', 'n_contigs', 'genome_mb'],
    'compatible-solute synthesis': ['ectA', 'ectB', 'ectC', 'otsA', 'otsB', 'betA', 'betB'],
    'compatible-solute uptake': ['betL', 'opuAA', 'opuAB', 'opuAC', 'opuBA', 'opuBB', 'opuBC',
                                 'opuCA', 'opuCB', 'opuCC', 'opu_permease', 'proP', 'proX', 'teaA'],
    'K+ uptake': ['trkA', 'trkH', 'ktrA', 'ktrB', 'ktrD', 'kdpA'],
    'Na+/H+ antiport': ['mrpA', 'mrpA_mrpB_fused', 'mrpB', 'mrpC', 'mrpD', 'mrpE', 'mrpF',
                        'mrpG', 'nhaA'],
    'mechanosensitive channels': ['mscL', 'mscS'],
    'other stress / housekeeping': ['cspA', 'galE', 'murB', 'mazG', 'gor', 'gshA', 'gshB', 'gshF'],
    'breadth summary': ['n_families_present', 'total_osmo_copies'],
}
EXAMPLES = {  # bacdive_id: label
    '1172': 'Bacillus subtilis',
    '12758': 'Pseudomonas aeruginosa',
    '6011': 'Halomonas elongata',
    '5998': 'Chromohalobacter salexigens',
    '133820': 'Halomonas ramblicola',
}


def group_of(f):
    for g, members in GROUPS.items():
        if f in members:
            return g
    return 'systems' if f.startswith('sys_') else 'other'


def shap_values(booster, X, names):
    d = xgb.DMatrix(X, feature_names=names)
    c = booster.predict(d, pred_contribs=True)
    return c[:, :-1], c[:, -1]


def direction(x, s):
    if np.nanstd(x) == 0 or np.nanstd(s) == 0:
        return np.nan
    return stats.spearmanr(x, s)[0]


def fmt(v, nd=3):
    return '' if v is None or (isinstance(v, float) and np.isnan(v)) else f'{v:.{nd}f}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--build', required=True)
    ap.add_argument('--pv-dir', default=None)
    ap.add_argument('--family-ranking', default=None,
                    help="analyze_extended.py's family_ranking.tsv; default: <pv-dir>/family_ranking.tsv")
    ap.add_argument('--top', type=int, default=20)
    args = ap.parse_args()
    out = os.path.join(args.build, 'explain')
    os.makedirs(out, exist_ok=True)

    schema = json.load(open(os.path.join(args.build, 'feature_schema.json')))
    names = schema['features']
    booster = xgb.Booster()
    booster.load_model(os.path.join(args.build, schema['model_file']))
    df = pd.read_csv(os.path.join(args.build, 'training_table.tsv'), sep='\t',
                     dtype={'bacdive_id': str})
    A = df[(df['strategy'] == 'salt_out') & df['nacl_max_pct'].notna()].reset_index(drop=True)
    X = A[names].to_numpy(float)
    y = A['nacl_max_pct'].to_numpy(float)
    sv, bias = shap_values(booster, X, names)
    pred = np.clip(sv.sum(1) + bias, 0, None)

    print('[explain] global importance', flush=True)
    # global importance --------------------------------------------------------
    imp = {k: booster.get_score(importance_type=k) for k in ('gain', 'cover', 'weight')}
    rows = []
    for j, n in enumerate(names):
        rows.append({'feature': n, 'group': group_of(n),
                     'mean_abs_shap': np.abs(sv[:, j]).mean(),
                     'mean_shap_present': sv[X[:, j] > 0, j].mean() if (X[:, j] > 0).any() else np.nan,
                     'mean_shap_absent': sv[X[:, j] == 0, j].mean() if (X[:, j] == 0).any() else np.nan,
                     'direction_rho': direction(X[:, j], sv[:, j]),
                     'prevalence': (X[:, j] > 0).mean(),
                     **{k: imp[k].get(n, 0.0) for k in imp}})
    ti = pd.DataFrame(rows).sort_values('mean_abs_shap', ascending=False).reset_index(drop=True)
    ti['rank'] = np.arange(1, len(ti) + 1)
    ti.to_csv(os.path.join(out, 'shap_importance.tsv'), sep='\t', index=False, float_format='%.4f')

    groups = sorted({group_of(n) for n in names})
    gsv = {g: sv[:, [j for j, n in enumerate(names) if group_of(n) == g]].sum(1) for g in groups}
    gi = pd.DataFrame([{'group': g, 'n_features': sum(group_of(n) == g for n in names),
                        'mean_abs_shap': np.abs(v).mean(), 'sd_shap': v.std(),
                        'min_shap': v.min(), 'max_shap': v.max()} for g, v in gsv.items()]
                      ).sort_values('mean_abs_shap', ascending=False)
    gi.to_csv(os.path.join(out, 'group_importance.tsv'), sep='\t', index=False, float_format='%.4f')

    print('[explain] stability over genus', flush=True)
    # stability over genus-blocked folds ------------------------------------------
    p = {k: v for k, v in schema['xgboost_params'].items() if k not in ('objective', 'tree_method')}
    stab = {}
    gk = GroupKFold(n_splits=5, shuffle=True, random_state=42)
    for i, (tr, te) in enumerate(gk.split(X, groups=A['genus'])):
        m = xgb.XGBRegressor(objective='reg:squarederror', tree_method='hist', random_state=42,
                             n_jobs=8, **p).fit(X[tr], y[tr])
        b = m.get_booster()
        b.feature_names = names
        s, _ = shap_values(b, X[te], names)
        stab[f'fold{i+1}'] = np.abs(s).mean(0)
    st = pd.DataFrame(stab, index=names)
    ranks = st.rank(ascending=False)
    st = st.join(ranks.add_suffix('_rank'))
    st['rank_min'], st['rank_max'] = ranks.min(axis=1), ranks.max(axis=1)
    st = st.loc[ti['feature']]
    st.index.name = 'feature'
    st.to_csv(os.path.join(out, 'importance_stability.tsv'), sep='\t', float_format='%.4f')
    rc = ranks.corr(method='spearman').to_numpy()
    fold_rank_rho = rc[np.triu_indices_from(rc, 1)]
    top10 = [set(ranks[c].nsmallest(10).index) for c in ranks.columns]
    top10_core = set.intersection(*top10)

    print('[explain] interactions', flush=True)
    # interactions --------------------------------------------------------------
    # interaction SHAP costs ~depth^2 x features more than plain SHAP; 500
    # strains are plenty to rank pairs
    sub = np.random.default_rng(0).choice(len(X), size=min(500, len(X)), replace=False)
    inter = booster.predict(xgb.DMatrix(X[sub], feature_names=names), pred_interactions=True)
    mi = np.abs(inter[:, :-1, :-1]).mean(0)
    ir = []
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            ir.append({'feature_a': names[a], 'feature_b': names[b], 'mean_abs_interaction': 2 * mi[a, b]})
    it = pd.DataFrame(ir).sort_values('mean_abs_interaction', ascending=False).head(30)
    it.to_csv(os.path.join(out, 'interactions.tsv'), sep='\t', index=False, float_format='%.4f')
    main_eff = np.abs(np.diagonal(inter[:, :-1, :-1], axis1=1, axis2=2)).mean(0).sum()
    inter_share = 1 - main_eff / (main_eff + 2 * mi[np.triu_indices_from(mi, 1)].sum())

    print('[explain] ridge reference', flush=True)
    # ridge reference ---------------------------------------------------------------
    Xl = np.log1p(np.clip(X, 0, None))
    sc = StandardScaler().fit(Xl)
    r = RidgeCV(alphas=np.logspace(-2, 4, 25)).fit(sc.transform(Xl), y)
    rr = pd.DataFrame({'feature': names, 'ridge_coef_pct_per_sd': r.coef_})
    fr_path = args.family_ranking or os.path.join(args.pv_dir or '', 'family_ranking.tsv')
    if os.path.exists(fr_path):
        fr = pd.read_csv(fr_path, sep='\t').rename(columns={'family': 'feature'})
        rr = rr.merge(fr[['feature', 'rho', 'rho_given_baselines', 'within_genus_r']], on='feature', how='left')
    rr = rr.merge(ti[['feature', 'mean_abs_shap', 'direction_rho', 'rank']], on='feature')
    rr.sort_values('rank').to_csv(os.path.join(out, 'ridge_coefficients.tsv'), sep='\t',
                                  index=False, float_format='%.4f')
    agree = rr.dropna(subset=['direction_rho'])
    sign_agree = (np.sign(agree['ridge_coef_pct_per_sd']) == np.sign(agree['direction_rho'])).mean()

    print('[explain] examples', flush=True)
    # examples -----------------------------------------------------------------------
    ex_ids = [b for b in EXAMPLES if b in set(A['bacdive_id'])]
    lo_i, hi_i = int(np.argmin(pred)), int(np.argmax(pred))
    for i, lab in [(lo_i, 'lowest prediction'), (hi_i, 'highest prediction')]:
        if A.loc[i, 'bacdive_id'] not in ex_ids:
            ex_ids.append(A.loc[i, 'bacdive_id'])
    exrows = []
    idx = {b: i for i, b in enumerate(A['bacdive_id'])}
    for b in ex_ids:
        i = idx[b]
        order = np.argsort(-np.abs(sv[i]))
        exrows.append({'bacdive_id': b, 'species': A.loc[i, 'species'], 'true': y[i],
                       'pred': pred[i], 'bias': float(bias[i]),
                       **{f'top{k+1}': f'{names[j]}={X[i, j]:g} ({sv[i, j]:+.2f})'
                          for k, j in enumerate(order[:6])}})
    ex = pd.DataFrame(exrows)
    ex.to_csv(os.path.join(out, 'examples.tsv'), sep='\t', index=False, float_format='%.2f')

    print('[explain] plots', flush=True)
    # plots --------------------------------------------------------------------------
    top = ti.head(args.top)
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(top) + 1))
    col = ['#2a6fdb' if d > 0 else '#d1495b' if d < 0 else '#888' for d in top['direction_rho'].fillna(0)]
    ax.barh(top['feature'][::-1], top['mean_abs_shap'][::-1], color=col[::-1])
    ax.set_xlabel('mean |SHAP| (% NaCl)')
    ax.set_title('Feature importance (blue: more copies -> higher NaCl max; red: lower)')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'importance_bar.png'), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 0.3 * len(gi) + 1))
    ax.barh(gi['group'][::-1], gi['mean_abs_shap'][::-1], color='#555')
    ax.set_xlabel('mean |summed SHAP| (% NaCl)')
    ax.set_title('Importance per functional group')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'importance_groups.png'), dpi=150)
    plt.close(fig)

    # beeswarm
    rng = np.random.default_rng(1)
    fig, ax = plt.subplots(figsize=(8, 0.35 * len(top) + 1))
    for k, f in enumerate(top['feature'][::-1]):
        j = names.index(f)
        x = X[:, j]
        cval = stats.rankdata(x) / len(x) if x.std() > 0 else np.full(len(x), 0.5)
        ax.scatter(sv[:, j], k + rng.uniform(-0.3, 0.3, len(x)), c=cval, cmap='coolwarm',
                   s=3, alpha=0.5, linewidths=0)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['feature'][::-1])
    ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel('SHAP value (% NaCl); colour: feature value rank (blue low, red high)')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'beeswarm.png'), dpi=150)
    plt.close(fig)

    # dependence plots for top 12
    t12 = list(ti['feature'].head(12))
    fig, axes = plt.subplots(3, 4, figsize=(14, 9))
    for ax, f in zip(axes.flat, t12):
        j = names.index(f)
        x = X[:, j]
        jit = rng.uniform(-0.15, 0.15, len(x)) if np.all(np.mod(x, 1) == 0) and x.max() < 40 else 0
        ax.scatter(x + jit, sv[:, j], s=3, alpha=0.3, color='#2a6fdb', linewidths=0)
        ax.axhline(0, color='k', lw=0.5)
        ax.set_title(f)
        ax.set_xlabel('value')
        ax.set_ylabel('SHAP (% NaCl)')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'dependence_top12.png'), dpi=130)
    plt.close(fig)

    # interaction heatmap of top 12
    ii = [names.index(f) for f in t12]
    M = 2 * mi[np.ix_(ii, ii)]
    np.fill_diagonal(M, np.nan)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(M, cmap='viridis')
    ax.set_xticks(range(len(t12)))
    ax.set_xticklabels(t12, rotation=90)
    ax.set_yticks(range(len(t12)))
    ax.set_yticklabels(t12)
    fig.colorbar(im, label='mean |SHAP interaction| (% NaCl)')
    ax.set_title('Pairwise interactions, top-12 features')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'interactions_top12.png'), dpi=150)
    plt.close(fig)

    # waterfalls
    ncol = 4
    nrow = int(np.ceil(len(ex_ids) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 4.2 * nrow), squeeze=False)
    for ax in axes.flat[len(ex_ids):]:
        ax.axis('off')
    for ax, b in zip(axes.flat, ex_ids):
        i = idx[b]
        order = np.argsort(-np.abs(sv[i]))[:8]
        rest = sv[i].sum() - sv[i, order].sum()
        labels = [f'{names[j]}={X[i, j]:g}' for j in order] + ['other features']
        vals = list(sv[i, order]) + [rest]
        start = float(bias[i])
        for k, (lab, v) in enumerate(zip(labels, vals)):
            ax.barh(len(vals) - k, v, left=start, color='#2a6fdb' if v > 0 else '#d1495b')
            start += v
        ax.set_yticks(range(len(vals), 0, -1))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(bias[i], color='k', lw=0.5, ls=':')
        ax.set_title(f'{A.loc[i, "species"]}\ntrue {y[i]:g}%, pred {pred[i]:.1f}%', fontsize=9)
        ax.set_xlabel('% NaCl')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'waterfall_examples.png'), dpi=130)
    plt.close(fig)

    print('[explain] README', flush=True)
    # README ------------------------------------------------------------------------
    cv = pd.read_csv(os.path.join(args.build, 'cv_results.tsv'), sep='\t')
    cvb = pd.read_csv(os.path.join(args.build, 'cv_by_bin.tsv'), sep='\t')
    g, f = schema['cv']['genus_blocked'], schema['cv']['family_blocked']
    bg, bf = schema['cv']['baseline_genus_blocked'], schema['cv']['baseline_family_blocked']
    bl_names = ' + '.join(f'`{x}`' for x in schema['cv']['baseline_features'])
    det = schema['salt_in_detector']
    dthr = [r for r in det['cv_genus_blocked'] if r['threshold'] == det['threshold']][0]
    sip = os.path.join(args.build, 'salt_in_predictions.tsv')
    si = pd.read_csv(sip, sep='\t') if os.path.exists(sip) else None

    def row(n):
        d = ti.set_index('feature').loc[n]
        arrow = 'higher' if d['direction_rho'] > 0.1 else 'lower' if d['direction_rho'] < -0.1 else 'mixed'
        return (f'| {int(d["rank"])} | `{n}` | {d["group"]} | {d["mean_abs_shap"]:.2f} | '
                f'{arrow} ({d["direction_rho"]:+.2f}) | {d["prevalence"]:.2f} | '
                f'{st.loc[n, "rank_min"]:.0f}–{st.loc[n, "rank_max"]:.0f} |')

    abl = cv[cv['experiment'].isin(['ablation', 'baseline', 'target=sqrt',
                                    'weights=1/sqrt(n_genus)', 'train without censored'])]
    abl_lines = '\n'.join(
        f'| {r.experiment} | {r.model} | {r.features} | {r.block} | {r.spearman:.3f} | {r.mae:.2f} | {r.r2:.3f} |'
        for r in abl.itertuples() if not pd.isna(r.spearman))
    bin_lines = '\n'.join(f'| {r.block} | {r.true_bin} | {r.n} | {r.mean_true:.1f} | {r.mean_pred:.1f} | '
                          f'{r.mae:.2f} | {r.bin_accuracy:.2f} |' for r in cvb.itertuples())
    iv_lines = '\n'.join(
        f'| {"" if b["pred_gt"] is None else b["pred_gt"]}–{"" if b["pred_le"] is None else b["pred_le"]} | '
        f'{b["n"]} | {b["resid_lo"]:+.2f} | {b["resid_hi"]:+.2f} |' for b in schema['intervals']['bins'])
    ex_lines = '\n'.join(f'| {r.species} ({r.bacdive_id}) | {r.true:g} | {r.pred:.1f} | {r.top1}; {r.top2}; {r.top3} |'
                         for r in ex.itertuples())
    gi_lines = '\n'.join(f'| {r.group} | {r.n_features} | {r.mean_abs_shap:.2f} |' for r in gi.itertuples())
    it_lines = '\n'.join(f'| `{r.feature_a}` × `{r.feature_b}` | {r.mean_abs_interaction:.3f} |'
                         for r in it.head(8).itertuples())

    fgroups = schema.get('family_groups') or {}
    not_rep = schema.get('families_not_reported') or {}
    if schema.get('annotate_report_levels'):
        fam_line = (f"the copy number of each of {len(schema['families'])} osmo families that "
                    "`osmotool annotate` reports per gene")
        fam_line += ''.join(f";\n  - `{g}`, the summed copy number of {'/'.join(ms)}, which annotate reports "
                            "only as this group" for g, ms in fgroups.items())
        fam_line += (f".\n\n  Families annotate does not report ({len(not_rep)}: "
                     f"{', '.join(sorted(not_rep))}) are left out: their labels can't be called reliably, "
                     "so the model doesn't read them")
    else:
        fam_line = f"the copy number of each of {len(schema['families'])} osmo families"
    card = f"""# osmo_refdb NaCl-max predictor ({schema['refdb_release']})

## What it predicts
The highest NaCl concentration (% w/v) at which a **salt-out** prokaryote grows,
as BacDive records it, from the output of `osmotool annotate`
(`--method {schema['annotate']['method']}`, Prodigal `{schema['annotate']['prodigal']}`) run with
osmo_refdb **{schema['refdb_release']}**. The model uses no taxonomy.

## How it works
- **Algorithm.** Gradient-boosted regression trees (XGBoost {schema['software']['xgboost']},
  `{schema['xgboost_params']['n_estimators']}` trees of depth ≤ {schema['xgboost_params']['max_depth']},
  learning rate {schema['xgboost_params']['learning_rate']:.3f}). The hyperparameters come from a
  random search scored by genus-blocked CV.
- **Combining the trees.** Each tree splits on gene copy counts and genome properties. The prediction is
  the sum of all trees: a weighted vote over thresholds such as "ectC ≥ 1" or "kdpA = 0".
- **Inputs.** {len(names)} features:
  - genome size and composition (`total_proteins`, `gc`, `n_contigs`);
  - {fam_line}.

  Breadth summaries, operon-level "systems" calls and `genome_mb` were tested and dropped: they add
  nothing over the family counts (see the ablations below).
- **Explanations.** Every prediction can be explained exactly. TreeSHAP splits each prediction into
  the training mean ({bias[0]:.2f}% NaCl) plus one contribution per feature, in % NaCl.
  `predict.py` and osmotool report the largest contributions per genome.
- **Training data.** {schema['training']['n_strains']:,} BacDive strains with a quantitative NaCl max and
  strategy salt-out, from {schema['training']['n_genera']:,} genera and {schema['training']['n_families']} families.
  The median label is {schema['training']['target_quantiles']['0.5']:g}% (IQR
  {schema['training']['target_quantiles']['0.25']:g}–{schema['training']['target_quantiles']['0.75']:g}%).

## How well it works (nested CV)
The test folds contain only genera (or families) the model never saw. The hyperparameter search
was rerun inside each training fold.

| model | blocking | Spearman | Pearson | MAE (% NaCl) | RMSE | R² |
|---|---|---|---|---|---|---|
| final ({len(names)} features) | genus | {g['spearman']:.3f} | {g['pearson']:.3f} | {g['mae']:.2f} | {g['rmse']:.2f} | {g['r2']:.3f} |
| final ({len(names)} features) | family | {f['spearman']:.3f} | {f['pearson']:.3f} | {f['mae']:.2f} | {f['rmse']:.2f} | {f['r2']:.3f} |
| baseline ({bl_names}) | genus | {bg['spearman']:.3f} | {bg['pearson']:.3f} | {bg['mae']:.2f} | {bg['rmse']:.2f} | {bg['r2']:.3f} |
| baseline ({bl_names}) | family | {bf['spearman']:.3f} | {bf['pearson']:.3f} | {bf['mae']:.2f} | {bf['rmse']:.2f} | {bf['r2']:.3f} |

The baseline model is the same algorithm and the same nested procedure, with its own
hyperparameter search, given only genome size and GC. The difference between the two is what the
osmo genes add:
- genus-blocked R² goes from {bg['r2']:.2f} to {g['r2']:.2f};
- family-blocked R² goes from {bf['r2']:.2f} to {f['r2']:.2f}.

For reference, predicting the global mean for every strain gives MAE
{cv[(cv.model == 'predict global mean') & (cv.block == 'genus')].mae.iloc[0]:.2f}.

Per true NaCl-max bin, from the nested OOF predictions:

| blocking | true bin | n | mean true | mean pred | MAE | bin accuracy |
|---|---|---|---|---|---|---|
{bin_lines}

The model **compresses the range**: it overestimates tolerant-poor strains and underestimates
extreme halotolerants. This is typical of regression to the mean when the features explain about
half of the variance. Read the prediction as a ranking and a rough level, not an exact value.

### Prediction interval
The interval is split-conformal: the 10th–90th percentile of the out-of-fold residuals, per
predicted-value bin.

| predicted % | n | residual q10 | residual q90 |
|---|---|---|---|
{iv_lines}

Coverage:
- {schema['intervals']['coverage']['genus_oof_in_sample']:.1%} on the genus-blocked OOF (calibrated on it);
- {schema['intervals']['coverage']['family_oof_with_genus_intervals']:.1%} on the family-blocked OOF;
- nominal 80%.

## What the model weighs
"Direction" is the Spearman correlation between a feature's value and its SHAP contribution:
"higher" means more copies push the prediction up. "Fold ranks" is the range of the feature's rank
across the five genus-blocked CV refits.

| rank | feature | group | mean \\|SHAP\\| (% NaCl) | direction | prevalence | fold ranks |
|---|---|---|---|---|---|---|
{chr(10).join(row(n) for n in ti['feature'].head(args.top))}

Summed per functional group:

| group | features | mean \\|SHAP\\| |
|---|---|---|
{gi_lines}

How stable these rankings are:
- Across the five genus-blocked refits, Spearman between the feature rankings is
  {fold_rank_rho.min():.2f}–{fold_rank_rho.max():.2f}.
- These {len(top10_core)} features are in the top 10 of every fold: {', '.join(f'`{x}`' for x in ti['feature'] if x in top10_core)}.

Interactions:
- They carry about {inter_share:.0%} of the total attribution; the rest is main effects.
- The largest pairs:

| pair | mean \\|interaction\\| |
|---|---|
{it_lines}

The linear reference (standardised ridge on log1p features) gives the same sign as the SHAP
direction for {sign_agree:.0%} of features. Its coefficients are in `explain/ridge_coefficients.tsv`,
next to the within-genus correlations from the BacDive validation (`family_ranking.tsv`).

Plots:
- `explain/importance_bar.png`, `explain/importance_groups.png` and `explain/beeswarm.png`: global importance.
- `explain/dependence_top12.png`: how each top feature's value maps to its contribution.
- `explain/interactions_top12.png`: interactions.
- `explain/waterfall_examples.png`: per-strain explanations.

### Example strains
These are training strains, so their predictions are in-sample.

| strain | BacDive % | predicted % | largest contributions (value, SHAP) |
|---|---|---|---|
{ex_lines}

## Ablations and sensitivity
Hyperparameters were fixed to the final ones, so read these as comparisons, not as honest estimates.

| experiment | model | features | blocking | Spearman | MAE | R² |
|---|---|---|---|---|---|---|
{abl_lines}

## Salt-in organisms (out of domain)
Salt-in organisms (Halobacteria, Salinibacter, Halanaerobiales) reach high NaCl tolerance by
accumulating KCl. They do not use the osmolyte genes this panel measures, so the model was not
trained on them.
"""
    if si is not None:
        card += (f"- **Out-of-domain readout.** On the {len(si)} salt-in BacDive strains the median true value "
                 f"is {si.nacl_max_pct.median():.1f}% and the median prediction {si.pred.median():.1f}%. "
                 f"The model underestimates them badly, as expected.\n")
    card += f"""- **Salt-in warning.** osmotool raises a warning from a separate logistic-regression detector
  (salt-in vs salt-out on the same features, {det['n_salt_in']} vs {det['n_salt_out']} strains)
  when p(salt-in) ≥ {det['threshold']}. In genus-blocked CV at that threshold it has precision
  {dthr['precision']:.2f} and recall {dthr['recall']:.2f} ({dthr['tp']} TP, {dthr['fp']} FP, {dthr['fn']} FN).
- **What the warning really flags.** Most false alarms are non-halophilic archaea and anaerobes
  that also carry few osmolyte genes. Read the warning as "this genome looks like the salt-in
  lineages; the NaCl prediction is not reliable", not as proof of a salt-in strategy.

## Limitations
- **Label quality.** BacDive's NaCl max is the highest concentration a lab tested with growth.
  - Values cluster on the test grid (0, 2, 5, 10, 15%).
  - They are capped at the highest concentration each lab tried, a cap that is not recorded.
  - {int(A['censored'].fillna(0).sum())} labels carry a "<" or ">" prefix. Dropping them from training changed nothing.
- **Taxonomy.** Genus alone explains ~75% of the label variance in BacDive. The model sees no
  taxonomy. Within a genus it captures only the part carried by osmo-gene content and genome size;
  across genera it effectively infers the lineage from gene content.
- **Tied to {schema['refdb_release']}.** Family counts depend on the refdb's HMM cutoffs; for
  example, ktrB prevalence fell from 0.43 to 0.08 between v9 and v10. Use the model only with
  `osmotool annotate --method hmm` against osmo_refdb {schema['refdb_release']}.
- **Assembly quality.** Fragmented or incomplete assemblies lose genes, which lowers the counts
  and therefore the prediction. `predict.py` flags genomes with more than
  {schema['applicability']['n_contigs_warn']:.0f} contigs (the 99th percentile in training), and genomes
  with more than {schema['applicability']['max_features_outside']} features outside the
  training 1st–99th percentile.
- **Associations, not mechanism.** SHAP shows what the model uses, not what causes salt tolerance.
  Correlated families (ectA/B/C, mrpA–G) share credit between them.
"""
    open(os.path.join(out, 'README.md'), 'w').write(card)
    print(ti.head(args.top)[['rank', 'feature', 'group', 'mean_abs_shap', 'direction_rho', 'prevalence']]
          .to_string(index=False, float_format='%.3f'))
    print(gi.to_string(index=False, float_format='%.3f'))
    print(f'fold-rank rho {fold_rank_rho.min():.2f}-{fold_rank_rho.max():.2f}; top10 core {sorted(top10_core)}')
    print(f'interaction share {inter_share:.2f}; ridge/SHAP sign agreement {sign_agree:.2f}')
    print(ex.to_string(index=False))
    print(f'wrote {out}/')


if __name__ == '__main__':
    main()
