#!/usr/bin/env python3
"""
12_package_markers.py -- copy a separately-built marker-gene release
(markers.yaml) into a main osmo_refdb release as its markers/ subdirectory,
for `osmotool profile`'s copies-per-genome normalisation.

Why a separate build at all: the universal single-copy markers must not sit
in the same DIAMOND db as the osmo families, or they would compete in
select_best_hits' best-alignment contest and shift the osmo families'
benchmarked precision/recall. So they are built and benchmarked on their own
(FAMILIES_FILE=markers.yaml bash run_pipeline.sh all <name>-markers) and
osmotool searches them in a second, independent DIAMOND pass over the same
reads. This step only packages; it never rebuilds anything.

Usage:
  python pipeline/12_package_markers.py \\
      --marker-release releases/v8-markers --release releases/v8

Writes releases/<release>/markers/:
  osmo_markers.dmnd                    DIAMOND db (marker TRAIN refs + decoys)
  osmo_markers.train_refs.faa          its source FASTA -- osmotool reads
                                       reference lengths from it (hit-weighted
                                       mean subject length per marker)
  osmo_markers.excluded_families.txt   decoy labels; never counted
  osmo_markers.manifest.tsv            one row per marker: n_train_refs,
                                       median/mean length (aa), plus the marker
                                       build's DIAMOND benchmark precision/
                                       recall/F1 when available
  osmo_markers.source.txt              provenance (source release, markers.yaml
                                       checksum, UniProt release if recorded)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import yaml

PREFIX = "osmo_markers"
SRC_PREFIX = "osmo_refdb"   # run_pipeline.sh's fixed output prefix


def fasta_lengths(path: Path) -> dict[str, list[int]]:
    """{family_label: [length_aa, ...]} from a family|accession|... FASTA."""
    lengths: dict[str, list[int]] = defaultdict(list)
    fam, n = None, 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if fam is not None:
                    lengths[fam].append(n)
                fam, n = line[1:].split("|", 1)[0].strip(), 0
            else:
                n += len(line.strip())
    if fam is not None:
        lengths[fam].append(n)
    return lengths


def diamond_metrics(summary: Path) -> dict[str, dict[str, str]]:
    """{family: {precision, recall, f1}} for method == diamond, if the marker
    build was benchmarked. Tolerant of a missing file or extra columns."""
    if not summary.is_file():
        return {}
    out: dict[str, dict[str, str]] = {}
    with summary.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("method", "").lower() != "diamond":
                continue
            fam = row.get("family")
            if fam:
                out[fam] = {k: row.get(k, "") for k in ("precision", "recall", "f1")}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--marker-release", required=True, type=Path,
                    help="release dir built from markers.yaml, e.g. releases/v8-markers")
    ap.add_argument("--release", required=True, type=Path,
                    help="main osmo_refdb release dir to package into, e.g. releases/v8")
    ap.add_argument("--markers", type=Path, default=Path("markers.yaml"))
    args = ap.parse_args()

    src, dst_release = args.marker_release, args.release
    dmnd = src / f"{SRC_PREFIX}.dmnd"
    faa = src / f"{SRC_PREFIX}.train_refs.faa"
    excluded = src / f"{SRC_PREFIX}.profile_excluded_families.txt"
    for p in (dmnd, faa, excluded):
        if not p.is_file():
            sys.exit(f"missing {p} -- run FAMILIES_FILE={args.markers} bash run_pipeline.sh build "
                     f"{src.name} first")
    if not (dst_release / f"{SRC_PREFIX}.dmnd").is_file():
        sys.exit(f"{dst_release} doesn't look like a main osmo_refdb release "
                 f"(no {SRC_PREFIX}.dmnd)")

    markers = [f["name"] for f in yaml.safe_load(args.markers.read_text())["families"]]
    decoys = [l.strip() for l in excluded.read_text().splitlines() if l.strip()]
    lengths = fasta_lengths(faa)

    missing = [m for m in markers if not lengths.get(m)]
    if missing:
        sys.exit(f"markers with no TRAIN references in {faa}: {', '.join(missing)}")
    unexpected = sorted(set(lengths) - set(markers) - set(decoys))
    if unexpected:
        sys.exit(f"{faa} contains labels that are neither markers nor excluded decoys: "
                 f"{', '.join(unexpected)}")

    metrics = diamond_metrics(src / "results" / "metrics" / "summary.tsv")
    if not metrics:
        print(f"WARNING: no DIAMOND benchmark metrics under {src}/results/metrics/ -- "
              f"run `benchmark` on the marker release before shipping it", file=sys.stderr)

    out = dst_release / "markers"
    out.mkdir(exist_ok=True)
    shutil.copy2(dmnd, out / f"{PREFIX}.dmnd")
    shutil.copy2(faa, out / f"{PREFIX}.train_refs.faa")
    (out / f"{PREFIX}.excluded_families.txt").write_text("".join(f"{d}\n" for d in decoys))

    with (out / f"{PREFIX}.manifest.tsv").open("w") as fh:
        fh.write("marker\tn_train_refs\tmedian_length_aa\tmean_length_aa\t"
                 "diamond_precision\tdiamond_recall\tdiamond_f1\n")
        for m in markers:
            ls = lengths[m]
            mt = metrics.get(m, {})
            fh.write(f"{m}\t{len(ls)}\t{statistics.median(ls):.0f}\t{statistics.mean(ls):.1f}\t"
                     f"{mt.get('precision', '')}\t{mt.get('recall', '')}\t{mt.get('f1', '')}\n")

    sha = hashlib.sha256(args.markers.read_bytes()).hexdigest()
    fetch_manifest = src / "refs" / "manifest.tsv"
    provenance = {}
    if fetch_manifest.is_file():
        with fetch_manifest.open() as fh:
            row = next(csv.DictReader(fh, delimiter="\t"), None) or {}
        provenance = {k: row[k] for k in ("uniprot_release", "date_fetched") if row.get(k)}
    with (out / f"{PREFIX}.source.txt").open("w") as fh:
        fh.write(f"marker_release\t{src.name}\n")
        fh.write(f"markers_yaml_sha256\t{sha}\n")
        for k, v in provenance.items():
            fh.write(f"{k}\t{v}\n")

    print(f"Packaged {len(markers)} markers (+{len(decoys)} decoy labels) into {out}/")


if __name__ == "__main__":
    main()
