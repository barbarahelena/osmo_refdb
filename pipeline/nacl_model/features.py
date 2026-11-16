#!/usr/bin/env python3
"""Genome features for the NaCl-max model, built from osmotool annotate output.

Single source of truth for what the model sees. Every feature comes from the
files `osmotool annotate` writes for one genome:

  <prefix>.gene_counts.tsv   per-family raw_count, plus the header lines
                             '# total_proteins' and (osmotool >= 0.8.0)
                             '# genome_bp', '# gc', '# n_contigs'
  <prefix>.systems.tsv       per-system present / subunit_count

No taxonomy and no BacDive data enter the features, so the model can be
applied to any genome osmotool annotates.

Feature groups (names as they appear in feature_schema.json):

  base       total_proteins, gc, n_contigs, genome_mb
  families   one column per osmo family: the raw_count (number of proteins
             whose best hit is that family). A family absent from
             gene_counts.tsv counts as 0 (osmotool <= 0.8 writes teaA only
             when it is non-zero).
  summary    n_families_present (families with raw_count > 0) and
             total_osmo_copies (sum of raw_count), both over the families
             listed in the schema -- the 'breadth' of
             pipeline/phenotype_validation/build_matrix.py.
  systems    sys_<name>_present (0/1) and sys_<name>_subunits.

features_from_annotate() is deliberately stdlib-only: it is the reference the
osmotool implementation is tested against.

Two definitions differ between the training data and osmotool, see
README.md: genome_bp in the BacDive training set counts A/C/G/T only
(run_random_panel.sh), osmotool's header counts every sequence character.
"""
from __future__ import annotations

import csv
import os

DECOY_SUFFIX = '_decoy'


def read_gene_counts(path):
    """-> ({header_key: value_str}, {family: raw_count})"""
    header, counts = {}, {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#'):
                parts = line[1:].strip().split('\t')
                if len(parts) >= 2:
                    header[parts[0].strip()] = parts[1].strip()
            elif line.startswith('gene\t') or not line:
                continue
            else:
                f = line.split('\t')
                if len(f) >= 2 and not f[0].endswith(DECOY_SUFFIX):
                    counts[f[0]] = float(f[1])
    return header, counts


def read_systems(path):
    """-> {system: (present 0/1, subunit_count, fragmented 0/1)}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r['system']] = (int(r['present'] == 'True'),
                                int(r['subunit_count']),
                                int(r.get('fragmented_assembly_flag') == 'True'))
    return out


def _num(x):
    if x is None or x == '' or x == 'NA':
        return None
    return float(x)


def features_from_annotate(prefix, families, systems, genome_stats=None):
    """Raw feature dict for one genome.

    prefix        annotate --out_prefix (files are <prefix>.gene_counts.tsv ...)
    families      the family columns to emit (from the schema)
    systems       the system names to emit (from the schema)
    genome_stats  optional {genome_bp, gc, n_contigs}; overrides the
                  gene_counts.tsv header (the BacDive training set predates
                  the osmotool 0.8.0 header and supplies them separately)

    Base features that cannot be found are None (the model treats them as
    missing). Also returns 'fragmented_systems', the number of systems with
    fragmented_assembly_flag set, for the applicability flags (not a model
    feature).
    """
    header, counts = read_gene_counts(f'{prefix}.gene_counts.tsv')
    sysd = read_systems(f'{prefix}.systems.tsv')
    gs = dict(genome_stats or {})
    bp = _num(gs.get('genome_bp', header.get('genome_bp')))
    feat = {
        'total_proteins': _num(header.get('total_proteins')),
        'gc': _num(gs.get('gc', header.get('gc'))),
        'n_contigs': _num(gs.get('n_contigs', header.get('n_contigs'))),
        'genome_mb': bp / 1e6 if bp is not None else None,
    }
    for f in families:
        feat[f] = counts.get(f, 0.0)
    feat['n_families_present'] = float(sum(1 for f in families if counts.get(f, 0) > 0))
    feat['total_osmo_copies'] = float(sum(counts.get(f, 0.0) for f in families))
    for s in systems:
        p, n, _ = sysd.get(s, (0, 0, 0))
        feat[f'sys_{s}_present'] = float(p)
        feat[f'sys_{s}_subunits'] = float(n)
    feat['fragmented_systems'] = sum(v[2] for v in sysd.values())
    return feat


def feature_names(families, systems, groups=('base', 'families', 'summary', 'systems')):
    out = []
    if 'base' in groups:
        out += ['total_proteins', 'gc', 'n_contigs', 'genome_mb']
    if 'families' in groups:
        out += list(families)
    if 'summary' in groups:
        out += ['n_families_present', 'total_osmo_copies']
    if 'systems' in groups:
        for s in systems:
            out += [f'sys_{s}_present', f'sys_{s}_subunits']
    return out


# ---------------------------------------------------------------------------
# Training table (BacDive). Needs pandas.
# ---------------------------------------------------------------------------

PHENO = ['bacdive_id', 'species', 'genus', 'family', 'domain', 'strategy',
         'nacl_max_pct', 'censored', 'nogrowth_at_pct', 'asm_level']


def discover_families_systems(annot_dir):
    """Union of family names (non-decoy) and system names over all strains."""
    fams, syss = set(), set()
    for bid in os.listdir(annot_dir):
        gc = os.path.join(annot_dir, bid, f'{bid}.gene_counts.tsv')
        if os.path.exists(gc):
            fams |= set(read_gene_counts(gc)[1])
            syss |= set(read_systems(os.path.join(annot_dir, bid, f'{bid}.systems.tsv')))
    return sorted(fams), sorted(syss)


def build_training_table(pv_dir, families=None, systems=None):
    """One row per BacDive strain in <pv_dir>/feature_matrix.tsv.

    Phenotype columns from the matrix (which already applied build_matrix.py's
    1 Mb assembly filter); features rebuilt from <pv_dir>/annotate/<bid>/ with
    features_from_annotate(), genome stats from
    <pv_dir>/random_panel/strains/<bid>/<bid>.genome_stats.tsv.
    Returns (DataFrame, families, systems).
    """
    import pandas as pd

    fm = pd.read_csv(os.path.join(pv_dir, 'feature_matrix.tsv'), sep='\t',
                     dtype={'bacdive_id': str})
    annot = os.path.join(pv_dir, 'annotate')
    if families is None or systems is None:
        f2, s2 = discover_families_systems(annot)
        families = families if families is not None else f2
        systems = systems if systems is not None else s2

    rows = []
    for _, p in fm[PHENO].iterrows():
        bid = p['bacdive_id']
        gs_path = os.path.join(pv_dir, 'random_panel', 'strains', bid, f'{bid}.genome_stats.tsv')
        gs = None
        if os.path.exists(gs_path):
            with open(gs_path) as fh:
                gs = next(csv.DictReader(fh, delimiter='\t'))
        feat = features_from_annotate(os.path.join(annot, bid, bid), families, systems, gs)
        rows.append({**p.to_dict(), **feat})
    df = pd.DataFrame(rows)
    return df, list(families), list(systems)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='Write the rebuilt training table (sanity check).')
    ap.add_argument('--pv-dir', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    df, fams, syss = build_training_table(a.pv_dir)
    df.to_csv(a.out, sep='\t', index=False)
    print(f'{len(df)} strains, {len(fams)} families, {len(syss)} systems -> {a.out}')
