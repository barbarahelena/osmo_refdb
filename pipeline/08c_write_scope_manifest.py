#!/usr/bin/env python3
"""
08c_write_scope_manifest.py — list families scoped out of reported output.

families.yaml can mark a family:
  scope: annotate_only        excluded from `osmotool profile`'s reported
                               output only (e.g. murB: a near-universal
                               housekeeping gene included specifically so
                               `osmotool annotate` can check co-occurrence
                               with galE/mazG on an assembled genome, but
                               meaningless in short-read community
                               profiling). Still visible in `annotate`.
  decoy_from_negatives: true  its "<name>_decoy" DIAMOND references (see
                               08a_build_decoy_refs.py) must never be
                               reported in EITHER mode -- they exist purely
                               to win select_best_hits' contest away from a
                               mislabeled call, not to be counted as a gene
                               family themselves.

Both kinds of family are still built into the combined DIAMOND/HMM
databases as normal -- they're harmless (or actively useful, for decoys)
to search against. This writes the flat, dependency-free lists `osmotool
profile`/`annotate --exclude_families` read to filter them out of
gene_counts.tsv.

Usage:
  python 08c_write_scope_manifest.py --families families.yaml \\
      --profile-out releases/<name>/osmo_refdb.profile_excluded_families.txt \\
      --annotate-out releases/<name>/osmo_refdb.annotate_excluded_families.txt
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--profile-out", type=Path, required=True)
    ap.add_argument("--annotate-out", type=Path, required=True)
    args = ap.parse_args()

    with open(args.families) as fh:
        data = yaml.safe_load(fh)

    annotate_only = [fam["name"] for fam in data["families"] if fam.get("scope") == "annotate_only"]
    decoy_labels = [f"{fam['name']}_decoy" for fam in data["families"] if fam.get("decoy_from_negatives")]

    profile_excluded = annotate_only + decoy_labels
    annotate_excluded = decoy_labels

    args.profile_out.write_text("\n".join(profile_excluded) + ("\n" if profile_excluded else ""))
    args.annotate_out.write_text("\n".join(annotate_excluded) + ("\n" if annotate_excluded else ""))

    print(f"profile-excluded: {len(profile_excluded)} families {profile_excluded}")
    print(f"annotate-excluded: {len(annotate_excluded)} families {annotate_excluded}")
    print(f"Written: {args.profile_out}")
    print(f"Written: {args.annotate_out}")


if __name__ == "__main__":
    main()
