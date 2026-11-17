#!/usr/bin/env python3
"""
verify_reproducible.py -- compare two builds of the same release and say
what is identical.

The test of a reproducibility claim is to build twice (same commit, same
image, same inputs) and compare, not to argue. Compares, for two release
directories:
  * every hmms/*.hmm (byte for byte, ignoring the DATE and COM lines HMMER
    writes), the GA cutoffs inside them;
  * DIAMOND cutoffs, cascade config and report levels (comment lines skipped);
  * the DIAMOND database (byte for byte; reported but only a warning if it
    differs, as the file may carry build metadata);
  * results/metrics/summary_by_report_unit.tsv and annotate/overall.tsv
    (every number), if both releases have been benchmarked.
Exit status 1 if anything in the first group differs.

Usage:
  python bin/verify_reproducible.py releases/v10-run1 releases/v10-run2
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def hmm_sha(path: Path) -> str:
    """HMM checksum without HMMER's DATE and COM lines (timestamp and the
    hmmbuild command line, which carries file paths): the two lines in which
    two builds of one alignment always differ."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for line in fh:
            if not line.startswith((b"DATE", b"COM")):
                h.update(line)
    return h.hexdigest()


def ga(path: Path) -> str:
    m = re.search(r"^GA\s+(.*)$", path.read_text(), re.M)
    return m.group(1) if m else ""


def table(path: Path) -> list[list[str]]:
    with open(path) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", type=Path)
    ap.add_argument("b", type=Path)
    args = ap.parse_args()
    a, b = args.a, args.b
    bad = 0

    names_a = {p.name for p in (a / "hmms").glob("*.hmm")}
    names_b = {p.name for p in (b / "hmms").glob("*.hmm")}
    only = sorted(names_a ^ names_b)
    diff = [n for n in sorted(names_a & names_b) if hmm_sha(a / "hmms" / n) != hmm_sha(b / "hmms" / n)]
    ga_diff = [n for n in diff if ga(a / "hmms" / n) != ga(b / "hmms" / n)]
    print(f"HMMs: {len(names_a & names_b)} in both, {len(diff)} differ "
          f"({len(ga_diff)} with a different GA), {len(only)} in one release only")
    for n in diff[:12]:
        print(f"   differs: {n}" + (f"  GA {ga(a / 'hmms' / n)} vs {ga(b / 'hmms' / n)}"
                                    if n in ga_diff else ""))
    for n in only:
        print(f"   only in one: {n}")
    bad += len(diff) + len(only)

    for name in ("osmo_refdb.diamond_cutoffs.tsv", "osmo_refdb.profile_cascade.tsv",
                 "osmo_refdb.profile_report_levels.tsv"):
        pa, pb = a / name, b / name
        if not (pa.exists() and pb.exists()):
            print(f"{name}: missing in one release")
            bad += 1
            continue
        same = table(pa) == table(pb)
        print(f"{name}: {'identical' if same else 'DIFFERENT'}")
        bad += 0 if same else 1

    da, db = a / "osmo_refdb.dmnd", b / "osmo_refdb.dmnd"
    if da.exists() and db.exists():
        print(f"DIAMOND database: {'byte-identical' if sha(da) == sha(db) else 'differs (warning only)'}")

    for name in ("results/metrics/summary_by_report_unit.tsv", "results/annotate/overall.tsv"):
        pa, pb = a / name, b / name
        if pa.exists() and pb.exists():
            ta, tb = table(pa), table(pb)
            n = sum(1 for x, y in zip(ta, tb) if x != y) + abs(len(ta) - len(tb))
            print(f"{name}: {'identical' if n == 0 else f'{n} rows differ'}")
            bad += 1 if n else 0

    print("\nREPRODUCED" if not bad else f"\nNOT identical ({bad} differences in the must-match group)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
