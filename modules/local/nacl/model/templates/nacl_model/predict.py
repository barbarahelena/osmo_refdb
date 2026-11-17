#!/usr/bin/env python3
"""Apply the NaCl-max model to osmotool annotate output (reference implementation).

This is the specification, and the test oracle, for the osmotool version:
the same schema, the same features, the same output columns.

  python bin/nacl_model/predict.py <model_dir> <annotate_prefix> [...] [-o out.tsv]

<model_dir> holds feature_schema.json and the model file: a release's
nacl_model/ directory, or train.py's build directory. Each prefix is an
`osmotool annotate --out_prefix` (so <prefix>.gene_counts.tsv exists). The
gene_counts.tsv header must carry '# gc' and '# n_contigs', which osmotool
>= 0.8.0 writes. For older annotate output, pass --genome-stats with a TSV of
prefix-basename, gc, n_contigs.

Output columns, one row per genome:
  genome               prefix basename
  nacl_max_pred        predicted NaCl max (% w/v), clipped at 0
  nacl_max_lo/_hi      80% prediction interval (split-conformal)
                       No tolerance class is reported: the model compresses
                       the range (true <=1% strains average a 3.1% prediction,
                       true >15% strains 13.2%), so a class read off the
                       prediction would be wrong for most low and extreme
                       strains. Use the value as a ranking plus interval.
  p_salt_in            salt-in detector probability
  flags                comma-separated: salt_in_like, out_of_range(<n>),
                       many_contigs (more contigs than 99% of the training
                       genomes; contig count is not linked to larger errors
                       in CV, so this marks the edge of the training range,
                       not a bad prediction), few_proteins,
                       missing:<feature>; 'ok' if none
  base_value           SHAP bias (mean training prediction)
  top_contributions    the 6 largest SHAP contributions, "feature=value:+x.xx"
  shap_<feature>       (with --all-shap) every contribution
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

import numpy as np
import xgboost as xgb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import features as F  # noqa: E402


def load(model_dir):
    # xgboost < 3.1 loads 3.x models without error but misreads their base
    # score, shifting every prediction by several % NaCl
    if tuple(int(''.join(c for c in x if c.isdigit()) or 0)
             for x in xgb.__version__.split('.')[:2]) < (3, 1):
        sys.exit(f'xgboost >= 3.1 needed (found {xgb.__version__}): older versions '
                 'misread this model')
    schema = json.load(open(os.path.join(model_dir, 'feature_schema.json')))
    booster = xgb.Booster()
    booster.load_model(os.path.join(model_dir, schema['model_file']))
    return schema, booster


def salt_in_probability(det, feat):
    z = det['intercept']
    for n, m, s, c in zip(det['features'], det['mean'], det['scale'], det['coef']):
        v = feat.get(n)
        x = math.log1p(max(v, 0.0)) if v is not None else m
        z += c * (x - m) / s
    return 1 / (1 + math.exp(-z))


def interval(schema, pred):
    for b in schema['intervals']['bins']:
        lo = -math.inf if b['pred_gt'] is None else b['pred_gt']
        hi = math.inf if b['pred_le'] is None else b['pred_le']
        if lo < pred <= hi:
            return max(pred + b['resid_lo'], 0.0), pred + b['resid_hi']
    b = schema['intervals']['bins'][0]
    return max(pred + b['resid_lo'], 0.0), pred + b['resid_hi']


def predict(schema, booster, feats, n_top=6):
    """feats: list of raw feature dicts (features_from_annotate). -> list of result dicts"""
    names = schema['features']
    X = np.array([[np.nan if f.get(n) is None else f[n] for n in names] for f in feats], float)
    d = xgb.DMatrix(X, feature_names=names)
    contrib = booster.predict(d, pred_contribs=True)
    raw = contrib.sum(1)
    if schema['target']['transform'] == 'sqrt':
        pred = np.clip(raw, 0, None) ** 2
    else:
        pred = np.clip(raw, schema['target'].get('clip_min', 0.0), None)
    app = schema['applicability']
    out = []
    for i, f in enumerate(feats):
        flags = []
        p_si = salt_in_probability(schema['salt_in_detector'], f)
        if p_si >= schema['salt_in_detector']['threshold']:
            flags.append('salt_in_like')
        n_out = sum(1 for n in names if f.get(n) is not None and
                    not (app['feature_q01'][n] <= f[n] <= app['feature_q99'][n]))
        if n_out > app['max_features_outside']:
            flags.append(f'out_of_range({n_out})')
        if (f.get('n_contigs') or 0) > app['n_contigs_warn']:
            flags.append('many_contigs')
        if (f.get('total_proteins') or 0) < app['min_total_proteins']:
            flags.append('few_proteins')
        flags += [f'missing:{n}' for n in names if f.get(n) is None]
        lo, hi = interval(schema, float(pred[i]))
        sv = contrib[i, :-1]
        order = np.argsort(-np.abs(sv))[:n_top]
        out.append({
            'nacl_max_pred': float(pred[i]), 'nacl_max_lo': lo, 'nacl_max_hi': hi,
            'p_salt_in': p_si, 'flags': ','.join(flags) or 'ok',
            'base_value': float(contrib[i, -1]),
            'top_contributions': ';'.join(f'{names[j]}={X[i, j]:g}:{sv[j]:+.2f}' for j in order),
            'shap': dict(zip(names, sv.tolist())),
        })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('model_dir')
    ap.add_argument('prefixes', nargs='+')
    ap.add_argument('-o', '--out', default='-')
    ap.add_argument('--genome-stats', help='TSV: genome, gc, n_contigs[, genome_bp]')
    ap.add_argument('--all-shap', action='store_true')
    args = ap.parse_args()
    schema, booster = load(args.model_dir)

    gstats = {}
    if args.genome_stats:
        with open(args.genome_stats) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                gstats[r['genome']] = r
    feats, gnames = [], []
    for p in args.prefixes:
        g = os.path.basename(p)
        feats.append(F.features_from_annotate(p, schema['families'], schema['systems'], gstats.get(g)))
        gnames.append(g)
    res = predict(schema, booster, feats)

    cols = ['genome', 'nacl_max_pred', 'nacl_max_lo', 'nacl_max_hi', 'p_salt_in',
            'flags', 'base_value', 'top_contributions']
    if args.all_shap:
        cols += [f'shap_{n}' for n in schema['features']]
    fh = sys.stdout if args.out == '-' else open(args.out, 'w', newline='')
    fh.write(f'# model nacl_max, osmo_refdb {schema["refdb_release"]}, '
             f'annotate --method {schema["annotate"]["method"]}\n')
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(cols)
    for g, r in zip(gnames, res):
        row = [g, f'{r["nacl_max_pred"]:.2f}', f'{r["nacl_max_lo"]:.2f}', f'{r["nacl_max_hi"]:.2f}',
               f'{r["p_salt_in"]:.3f}', r['flags'], f'{r["base_value"]:.3f}',
               r['top_contributions']]
        if args.all_shap:
            row += [f'{r["shap"][n]:.4f}' for n in schema['features']]
        w.writerow(row)
    if fh is not sys.stdout:
        fh.close()


if __name__ == '__main__':
    main()
