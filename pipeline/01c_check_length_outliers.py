#!/usr/bin/env python3
"""
01c_check_length_outliers.py — flag/drop sequences whose length is way off
from their family's typical length (likely multi-domain fusion proteins or
partial fragments, not clean single-domain family members).

Motivating case: benchmarking proX (v2 release) showed DIAMOND false
positives traced back to two "proX" positives (A0A3N6Q0R4, A0ABZ2UTQ7) that
are actually unreviewed, low-confidence automatic annotations at 641/659 aa
-- roughly double E. coli ProX's ~330 aa -- because they fuse the PF04069
substrate-binding domain to an unrelated PF00528 permease domain. The same
pattern showed up on the negative side: an ectC hard-negative (A0A918UXX5,
867 aa) turned out to carry a real PF04069 domain fused alongside its
Cupin_2 domain, so a read landing on that fused region legitimately (and
correctly) scored as proX-like -- not a DIAMOND error, a reference-set
contamination issue.

This check runs on both the positive and (already purity-filtered)
negative set per family: any sequence whose length falls outside
[median / max-ratio, median * max-ratio] is flagged as a likely fusion
protein (too long) or fragment (too short) and moved out of the set the
rest of the pipeline consumes, rather than silently discarded -- flagged
sequences are kept in a separate file for manual review.

Idempotent: the first run snapshots each as-is FASTA (refs/<family>.
<positive|negative>.faa) to a "pre_length_filter" copy and always re-filters
from that snapshot, so re-running with a different --max-ratio is safe.

Usage:
  python 01c_check_length_outliers.py --refs refs --families families.yaml \\
      --max-ratio 1.5
Output:
  refs/<family>.<label>.faa                    (filtered in place)
  refs/<family>.<label>.pre_length_filter.faa  (untouched snapshot)
  refs/<family>.<label>.length_outliers.faa    (excluded, for manual review)
  refs/length_outlier_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import statistics
from pathlib import Path

import yaml


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header, chunks = None, []
    if not path.exists():
        return records
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--max-ratio", type=float, default=1.5,
                     help="Flag sequences longer than median*max-ratio or "
                          "shorter than median/max-ratio")
    ap.add_argument("--min-sequences", type=int, default=5,
                     help="Skip a family/label if it has fewer sequences "
                          "than this -- too few to establish a meaningful "
                          "median length")
    args = ap.parse_args()

    families = load_family_names(args.families)
    manifest_rows = []
    # Written fresh, unconditionally, whenever upstream actually reran (01
    # for the fetch itself, 01b for purity-filtered negatives) -- unlike
    # faa_path (which this script overwrites every run), these are stable
    # signals for "has anything upstream changed since our last snapshot",
    # used below to detect a stale snapshot from before a families.yaml
    # edit. Using the later of the two for both labels is an intentional
    # simplification: occasionally re-snapshotting when strictly
    # unnecessary is harmless, silently using a stale snapshot is not.
    upstream_manifests = [args.refs / "manifest.tsv", args.refs / "negative_purity_manifest.tsv"]
    upstream_mtime = max((p.stat().st_mtime for p in upstream_manifests if p.exists()), default=None)

    for family in families:
        for label in ("positive", "negative"):
            faa_path = args.refs / f"{family}.{label}.faa"
            snapshot_path = args.refs / f"{family}.{label}.pre_length_filter.faa"
            outliers_path = args.refs / f"{family}.{label}.length_outliers.faa"

            if not faa_path.exists() or faa_path.stat().st_size == 0:
                print(f"[{family}/{label}] SKIP: no sequences fetched")
                continue

            # Always re-filter from the untouched snapshot, so this step is
            # safe to re-run (e.g. with a different --max-ratio) -- EXCEPT
            # if upstream (fetch or purity-filtering) genuinely reran more
            # recently than the snapshot, in which case the snapshot is
            # stale and must not silently shadow fresh data in faa_path.
            stale = (
                snapshot_path.exists()
                and upstream_mtime is not None
                and upstream_mtime > snapshot_path.stat().st_mtime
            )
            if stale:
                print(f"[{family}/{label}] Upstream fetch/purity data is newer than "
                      f"the existing snapshot -- re-snapshotting")
            if not snapshot_path.exists() or stale:
                shutil.copy(faa_path, snapshot_path)
            records = parse_fasta(snapshot_path)

            if len(records) < args.min_sequences:
                print(f"[{family}/{label}] SKIP: only {len(records)} sequences, "
                      f"too few for a meaningful median length")
                continue

            lengths = [len(seq) for _, seq in records]
            median_len = statistics.median(lengths)
            lower_bound = median_len / args.max_ratio
            upper_bound = median_len * args.max_ratio

            kept, flagged = [], []
            for header, seq in records:
                (kept if lower_bound <= len(seq) <= upper_bound else flagged).append((header, seq))

            write_fasta(faa_path, kept)
            write_fasta(outliers_path, flagged)

            print(f"[{family}/{label}] median={median_len:.0f}aa, "
                  f"allowed=[{lower_bound:.0f}, {upper_bound:.0f}]aa -> "
                  f"{len(kept)} kept, {len(flagged)} flagged")
            if flagged:
                worst = max(flagged, key=lambda r: abs(len(r[1]) - median_len))
                print(f"    e.g. {worst[0]} at {len(worst[1])}aa (median {median_len:.0f}aa) "
                      f"-- likely a fusion protein or fragment, see {outliers_path}")

            manifest_rows.append({
                "family": family,
                "label": label,
                "n_total": len(records),
                "n_kept": len(kept),
                "n_flagged": len(flagged),
                "median_length_aa": round(median_len, 1),
                "max_ratio": args.max_ratio,
            })

    if not manifest_rows:
        print("No families checked (no positive/negative ref FASTAs found).")
        return

    manifest_path = args.refs / "length_outlier_manifest.tsv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    n_with_flags = sum(1 for r in manifest_rows if r["n_flagged"] > 0)
    print(f"\nDone. {len(manifest_rows)} family/label sets checked, "
          f"{n_with_flags} had >=1 length outlier flagged.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()