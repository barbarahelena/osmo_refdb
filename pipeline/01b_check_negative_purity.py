#!/usr/bin/env python3
"""
01b_check_negative_purity.py — flag/drop hard negatives that are secretly
true positives.

Motivating bug (v1 release): betL's negative_query excluded only the
literal string "gene:betL", so BCCT-family orthologs discovered under a
different gene symbol in another organism (betP, opuD, betS -- all
genuine glycine-betaine transporters) sat in the "negative" set unnoticed.
Both DIAMOND and HMM collapsed to ~3% precision on betL as a result --
not because either method is bad, but because the ground truth was
self-contradictory. families.yaml has since been hand-fixed for betL, but
nothing stopped the same mistake from recurring silently for any other
family (existing or new), since gene-symbol string exclusion can never
enumerate every historical synonym a curator doesn't already know about.

This step makes that check automatic and family-agnostic: for every
family, align each hard-negative sequence against that family's full
positive set with DIAMOND and flag any negative that's suspiciously
similar (default: >=70% identity) as a likely mislabeled ortholog rather
than a true hard negative. Flagged sequences are removed from
refs/<family>.negative.faa (the file every downstream step reads) and
saved separately for manual review, not silently discarded.

Idempotent: the very first run snapshots the as-fetched negative set to
refs/<family>.negative.raw.faa and always re-filters from that snapshot,
so re-running this step (e.g. after lowering --threshold) never operates
on an already-filtered file.

Usage:
  python 01b_check_negative_purity.py --refs refs --families families.yaml \\
      --threshold 70.0
Output:
  refs/<family>.negative.faa          (filtered -- what downstream steps use)
  refs/<family>.negative.raw.faa      (untouched snapshot of the original fetch)
  refs/<family>.negative.flagged.faa  (excluded sequences, for manual review)
  refs/negative_purity_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path

import yaml


def load_family_names(path: Path) -> list[str]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [fam["name"] for fam in data["families"]]


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header, seq = None, []
    if not path.exists():
        return records
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def max_identity_per_negative(
    pos_path: Path, neg_path: Path, threads: int, tmp_dir: Path
) -> dict[str, float]:
    """Return {negative_header: best % identity to any positive sequence}.
    Negatives with no DIAMOND hit at all are absent from the dict (identity 0)."""
    db_path = tmp_dir / "positives.dmnd"
    hits_path = tmp_dir / "hits.tsv"

    subprocess.run(
        ["diamond", "makedb", "--in", str(pos_path), "--db", str(db_path), "--quiet"],
        check=True,
    )
    subprocess.run(
        ["diamond", "blastp",
         "--query", str(neg_path), "--db", str(db_path),
         "--out", str(hits_path), "--outfmt", "6", "qseqid", "pident",
         "--max-target-seqs", "1", "--threads", str(threads), "--quiet"],
        check=True,
    )

    best_identity: dict[str, float] = {}
    with open(hits_path) as fh:
        for line in fh:
            qseqid, pident = line.rstrip("\n").split("\t")
            pident = float(pident)
            if pident > best_identity.get(qseqid, -1.0):
                best_identity[qseqid] = pident
    return best_identity


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--threshold", type=float, default=70.0,
                     help="Max %% identity to any positive above which a "
                          "negative is flagged as a likely mislabeled ortholog")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    families = load_family_names(args.families)
    manifest_rows = []
    # 01_fetch_refs.py writes this fresh, unconditionally, every time it
    # runs -- unlike neg_path (which this script itself overwrites below),
    # it's a stable signal for "when did the fetch last genuinely run",
    # used to detect a stale raw snapshot from before a families.yaml
    # query edit (see snapshot staleness check below).
    fetch_manifest_path = args.refs / "manifest.tsv"

    for family in families:
        pos_path = args.refs / f"{family}.positive.faa"
        neg_path = args.refs / f"{family}.negative.faa"
        raw_path = args.refs / f"{family}.negative.raw.faa"
        flagged_path = args.refs / f"{family}.negative.flagged.faa"

        if not pos_path.exists() or pos_path.stat().st_size == 0:
            print(f"[{family}] SKIP: no positive set to compare against")
            continue
        if not neg_path.exists() or neg_path.stat().st_size == 0:
            print(f"[{family}] SKIP: no negative set fetched")
            continue

        # Always re-filter from the untouched as-fetched snapshot, so this
        # step is safe to re-run (e.g. with a different --threshold) --
        # EXCEPT if a fetch has genuinely happened more recently than the
        # snapshot (e.g. a negative_query edit in families.yaml), in which
        # case the old snapshot is stale and must not silently shadow the
        # fresh fetch sitting in neg_path. Comparing against neg_path's own
        # mtime wouldn't work here since this script overwrites neg_path
        # itself every run -- fetch_manifest_path is the stable reference.
        stale = (
            raw_path.exists()
            and fetch_manifest_path.exists()
            and fetch_manifest_path.stat().st_mtime > raw_path.stat().st_mtime
        )
        if stale:
            print(f"[{family}] Fetch is newer than the existing raw snapshot -- "
                  f"re-snapshotting (a families.yaml query edit since the last "
                  f"run would otherwise be silently ignored)")
        if not raw_path.exists() or stale:
            shutil.copy(neg_path, raw_path)
        source_records = parse_fasta(raw_path)

        with tempfile.TemporaryDirectory() as tmp:
            identities = max_identity_per_negative(pos_path, raw_path, args.threads, Path(tmp))

        kept, flagged = [], []
        for header, seq in source_records:
            pident = identities.get(header, 0.0)
            (flagged if pident >= args.threshold else kept).append((header, seq, pident))

        write_fasta(neg_path, [(h, s) for h, s, _ in kept])
        write_fasta(flagged_path, [(h, s) for h, s, _ in flagged])

        print(f"[{family}] {len(source_records)} negatives -> "
              f"{len(kept)} kept, {len(flagged)} flagged (>= {args.threshold}% identity to a positive)")
        if flagged:
            worst = max(flagged, key=lambda r: r[2])
            print(f"    e.g. {worst[0]} at {worst[2]:.1f}% identity to a positive -- "
                  f"likely a mislabeled ortholog, see {flagged_path}")

        manifest_rows.append({
            "family": family,
            "n_total_negatives": len(source_records),
            "n_kept": len(kept),
            "n_flagged": len(flagged),
            "threshold_pident": args.threshold,
        })

    if not manifest_rows:
        print("No families checked (no positive/negative ref pairs found).")
        return

    manifest_path = args.refs / "negative_purity_manifest.tsv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    n_with_flags = sum(1 for r in manifest_rows if r["n_flagged"] > 0)
    print(f"\nDone. {len(manifest_rows)} families checked, "
          f"{n_with_flags} had >=1 flagged negative.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
