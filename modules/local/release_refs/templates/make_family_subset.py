#!/usr/bin/env python3
"""
make_family_subset.py — write a families.yaml containing only the named
families, for a fast test-drive of one or two families through the whole
pipeline (fetch -> QC -> calibrate -> benchmark, ~minutes) before
committing to a full rebuild of every family (~90 minutes).

Generates the subset from the real families.yaml rather than having you
hand-maintain a separate copy, so it can't silently drift out of sync
(e.g. after editing a query in the real file and forgetting to mirror the
change).

Usage:
  python make_family_subset.py --families families.yaml \\
      --only galE,mazG,murB --out families_test.yaml

Then run the whole pipeline against just that subset, in a throwaway
release directory so the real release is never touched:

  FAMILIES_FILE=families_test.yaml bash run_pipeline.sh all v3-test
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--only", required=True,
                     help="Comma-separated family names to include, e.g. galE,mazG,murB")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    wanted = [name.strip() for name in args.only.split(",") if name.strip()]

    with open(args.families) as fh:
        data = yaml.safe_load(fh)

    all_names = {fam["name"] for fam in data["families"]}
    missing = [name for name in wanted if name not in all_names]
    if missing:
        raise SystemExit(f"Not found in {args.families}: {missing}")

    subset = {"families": [fam for fam in data["families"] if fam["name"] in wanted]}

    with open(args.out, "w") as fh:
        yaml.safe_dump(subset, fh, default_flow_style=False, sort_keys=False)

    print(f"{len(subset['families'])} families written to {args.out}: {wanted}")
    print(f"Run: FAMILIES_FILE={args.out} bash run_pipeline.sh all <test_release_name>")


if __name__ == "__main__":
    main()
