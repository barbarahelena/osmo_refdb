#!/usr/bin/env python3
"""
13_package_nacl_model.py -- copy a trained NaCl-max model (bin/nacl_model/
train.py + explain.py output) into a release as its nacl_model/ subdirectory,
for osmotool to apply to annotated genomes. Packaging only; never retrains.

The model is tied to the release whose annotate output it was trained on:
family counts depend on that release's HMM cutoffs. The script refuses a
model whose schema names another release.

Usage:
  python bin/13_package_nacl_model.py \\
      --build releases/v10/nacl_model_build --release releases/v10

Writes releases/<release>/nacl_model/:
  nacl_max.xgb.json        XGBoost model (JSON, load with xgboost.Booster)
  feature_schema.json      feature order and definitions, target transform,
                           prediction-interval table, applicability ranges,
                           salt-in detector coefficients, CV metrics,
                           provenance
  README.md                how the model works and what it weighs
  explain/*.png, *.tsv     the importance plots and tables behind the card
  predict.py, features.py  reference inference (osmotool's test oracle)
  self_check.json          5 training genomes' features and the model's raw
                           output for them, computed here; osmotool
                           recomputes them on load and refuses the model if
                           they differ (older xgboost versions load the
                           model without error but misread its base score)
  SHA256SUMS               checksums of the files above
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPLAIN_FILES = ['shap_importance.tsv', 'group_importance.tsv', 'importance_stability.tsv',
                 'interactions.tsv', 'ridge_coefficients.tsv', 'examples.tsv',
                 'importance_bar.png', 'importance_groups.png', 'beeswarm.png',
                 'dependence_top12.png', 'interactions_top12.png', 'waterfall_examples.png']
BUILD_FILES = ['cv_results.tsv', 'cv_by_bin.tsv', 'salt_in_predictions.tsv']


def write_self_check(build: Path, schema: dict, model_p: Path, out: Path, n: int = 5) -> None:
    """Raw model output for the first n training genomes, for osmotool to verify."""
    import csv
    import xgboost as xgb
    names = schema['features']
    with open(build / 'training_table.tsv') as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t')
                if r['strategy'] == 'salt_out' and r['nacl_max_pct']][:n]
    X = [[float(r[f]) for f in names] for r in rows]
    b = xgb.Booster()
    b.load_model(str(model_p))
    margin = b.predict(xgb.DMatrix(X, feature_names=names), output_margin=True)
    out.write_text(json.dumps({
        'description': 'raw model output (output_margin) for these feature rows; '
                       'a loader must reproduce it within tolerance',
        'xgboost_version': xgb.__version__, 'tolerance': 1e-4,
        'features': names, 'rows': X, 'output_margin': [float(v) for v in margin],
    }, indent=1))


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--build', required=True, type=Path, help="train.py's --out directory")
    ap.add_argument('--release', required=True, type=Path)
    ap.add_argument('--release-name', default=None, help="default: the --release directory's name")
    ap.add_argument('--out', type=Path, default=None,
                    help='write the packaged model here instead of <release>/nacl_model '
                         '(the Nextflow pipeline publishes it into the release)')
    args = ap.parse_args()

    schema_p = args.build / 'feature_schema.json'
    if not schema_p.is_file():
        sys.exit(f'missing {schema_p} -- run bin/nacl_model/train.py first')
    schema = json.loads(schema_p.read_text())
    model_p = args.build / schema['model_file']
    card_p = args.build / 'explain' / 'README.md'
    for p in (model_p, card_p):
        if not p.is_file():
            sys.exit(f'missing {p} -- run train.py and explain.py first')
    release_name = args.release_name or args.release.name
    if schema['refdb_release'] != release_name:
        sys.exit(f"model was trained on {schema['refdb_release']} annotate output, "
                 f'refusing to package it into {release_name}')
    if not (args.release / 'osmo_refdb.dmnd').is_file():
        sys.exit(f'{args.release} is not a built release (no osmo_refdb.dmnd)')

    dst = args.out or args.release / 'nacl_model'
    if dst.exists():
        shutil.rmtree(dst)
    (dst / 'explain').mkdir(parents=True)
    shutil.copy2(model_p, dst / model_p.name)
    shutil.copy2(schema_p, dst / schema_p.name)
    shutil.copy2(card_p, dst / 'README.md')
    for f in EXPLAIN_FILES:
        if (args.build / 'explain' / f).is_file():
            shutil.copy2(args.build / 'explain' / f, dst / 'explain' / f)
    for f in BUILD_FILES:
        if (args.build / f).is_file():
            shutil.copy2(args.build / f, dst / 'explain' / f)
    for f in ('predict.py', 'features.py'):
        shutil.copy2(HERE / 'nacl_model' / f, dst / f)
    write_self_check(args.build, schema, model_p, dst / 'self_check.json')

    files = sorted(p for p in dst.rglob('*') if p.is_file())
    (dst / 'SHA256SUMS').write_text(
        ''.join(f'{sha256(p)}  {p.relative_to(dst)}\n' for p in files))
    g = schema['cv']['genus_blocked']
    print(f'packaged {len(files)} files into {dst}')
    print(f"  {len(schema['features'])} features, n_train {schema['training']['n_strains']}, "
          f"genus-blocked Spearman {g['spearman']:.3f}, MAE {g['mae']:.2f}")


if __name__ == '__main__':
    main()
