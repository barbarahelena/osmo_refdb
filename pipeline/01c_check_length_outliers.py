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

Fusion-partner families (families.yaml: `fusion_partner: <other family>`,
e.g. mrpA/mrpD, whose gene product is a single fused ORF in some lineages
-- see families.yaml): a sequence roughly double the family's own median
length would normally get dropped here as a probable fusion artifact, which
is exactly wrong for these two -- the fused ORF is a real, valid detection
target (Task 1b), just one the standard single-family MSA/HMM can't
represent (aligning a ~1300aa fused sequence alongside ~800aa standalone
mrpA sequences would badly gap the alignment for everyone). So for these
families specifically, an over-length candidate is checked against a SECOND
window centered on (own_median + partner_median) instead of being flagged
as an outlier outright; if it also carries `fusion_marker_pfam` (Pfam
domain evidence fetched by 01_fetch_refs.py, in
refs/<family>.positive.domains.tsv, when available) it's routed to
refs/<family>.positive.fusion_candidates.faa instead of length_outliers.faa
-- excluded from the normal alignment/HMM path (like a true outlier) but
picked up separately by 08d_build_fusion_refs.py as a real, reportable
DIAMOND-only reference. Domain evidence is confirmatory, not required: if
domains.tsv isn't present (e.g. after re-running 01c alone without rerunning
01), the length-window match alone is enough to classify a candidate as a
fusion instead of discarding it -- a false positive here just means a
handful of oversized sequences reach 08d instead of the outlier pile, not a
lost detection target.

Usage:
  python 01c_check_length_outliers.py --refs refs --families families.yaml \\
      --max-ratio 1.5
Output:
  refs/<family>.<label>.faa                      (filtered in place)
  refs/<family>.<label>.pre_length_filter.faa    (untouched snapshot)
  refs/<family>.<label>.length_outliers.faa      (excluded, for manual review)
  refs/<family>.positive.fusion_candidates.faa   (fusion-partner families only)
  refs/length_outlier_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import statistics
from pathlib import Path

import yaml


def load_families(path: Path) -> list[dict]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return data["families"]


def load_family_names(path: Path) -> list[str]:
    return [fam["name"] for fam in load_families(path)]


def load_domain_evidence(refs_dir: Path, family: str) -> dict[str, list[str]]:
    """accession -> Pfam domain list, from 01_fetch_refs.py's fusion-partner
    domain fetch. Returns {} if not present (fusion detection then falls
    back to length evidence alone -- see module docstring)."""
    path = refs_dir / f"{family}.positive.domains.tsv"
    if not path.exists():
        return {}
    domains = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            domains[row["accession"]] = row["pfam_domains"].split(";") if row["pfam_domains"] else []
    return domains


def header_accession(header: str) -> str:
    """<tag>|UniProtID|Organism -> UniProtID (see 01_fetch_refs.py:retag_fasta)."""
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header


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

    fam_defs = load_families(args.families)
    families = [fam["name"] for fam in fam_defs]
    fusion_partner = {fam["name"]: fam["fusion_partner"] for fam in fam_defs if fam.get("fusion_partner")}
    fusion_marker = {fam["name"]: fam["fusion_marker_pfam"] for fam in fam_defs if fam.get("fusion_marker_pfam")}
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

    # Pass 1: snapshot (if needed) and compute each family/label's own
    # median length up front -- fusion-partner families need to know their
    # PARTNER's median before they can build the fusion-length window
    # below, and the partner may sort later in families.yaml.
    all_records: dict[tuple[str, str], list[tuple[str, str]]] = {}
    medians: dict[tuple[str, str], float] = {}
    for family in families:
        for label in ("positive", "negative"):
            faa_path = args.refs / f"{family}.{label}.faa"
            snapshot_path = args.refs / f"{family}.{label}.pre_length_filter.faa"

            if not faa_path.exists() or faa_path.stat().st_size == 0:
                continue

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
                continue
            all_records[(family, label)] = records
            medians[(family, label)] = statistics.median(len(seq) for _, seq in records)

    # Pass 2: filter, with fusion-length accommodation for declared pairs.
    for family in families:
        for label in ("positive", "negative"):
            faa_path = args.refs / f"{family}.{label}.faa"
            outliers_path = args.refs / f"{family}.{label}.length_outliers.faa"
            fusion_path = args.refs / f"{family}.positive.fusion_candidates.faa"

            if (family, label) not in all_records:
                if faa_path.exists() and faa_path.stat().st_size > 0:
                    print(f"[{family}/{label}] SKIP: only {len(parse_fasta(faa_path))} sequences, "
                          f"too few for a meaningful median length")
                else:
                    print(f"[{family}/{label}] SKIP: no sequences fetched")
                continue

            records = all_records[(family, label)]
            median_len = medians[(family, label)]
            lower_bound = median_len / args.max_ratio
            upper_bound = median_len * args.max_ratio

            # Fusion window: only meaningful for a positive set whose family
            # declares a fusion_partner AND that partner's own median is
            # known (from this same pass) -- see module docstring.
            partner = fusion_partner.get(family) if label == "positive" else None
            partner_median = medians.get((partner, "positive")) if partner else None
            fusion_lower = fusion_upper = None
            if partner_median is not None:
                fusion_center = median_len + partner_median
                fusion_lower = fusion_center / args.max_ratio
                fusion_upper = fusion_center * args.max_ratio
            marker_pfam = fusion_marker.get(family)
            domain_evidence = load_domain_evidence(args.refs, family) if partner else {}

            kept, flagged, fusion_hits = [], [], []
            for header, seq in records:
                n = len(seq)
                if lower_bound <= n <= upper_bound:
                    kept.append((header, seq))
                elif fusion_lower is not None and fusion_lower <= n <= fusion_upper:
                    accession = header_accession(header)
                    has_marker = (
                        not domain_evidence or not marker_pfam
                        or marker_pfam in domain_evidence.get(accession, [])
                    )
                    (fusion_hits if has_marker else flagged).append((header, seq))
                else:
                    flagged.append((header, seq))

            write_fasta(faa_path, kept)
            write_fasta(outliers_path, flagged)
            if partner:
                write_fasta(fusion_path, fusion_hits)

            print(f"[{family}/{label}] median={median_len:.0f}aa, "
                  f"allowed=[{lower_bound:.0f}, {upper_bound:.0f}]aa -> "
                  f"{len(kept)} kept, {len(flagged)} flagged"
                  + (f", {len(fusion_hits)} fusion candidates -> {fusion_path}" if partner else ""))
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
                "n_fusion_candidates": len(fusion_hits) if partner else "",
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