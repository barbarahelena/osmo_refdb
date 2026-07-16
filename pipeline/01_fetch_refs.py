#!/usr/bin/env python3
"""
01_fetch_refs.py — osmo_refdb reference sequence fetcher

Fetches, per gene family defined in families.yaml:
  * a broad POSITIVE set from UniProt (reviewed + unreviewed, all bacterial
    taxa) — used to build the MSA / HMM / DIAMOND db and to later hold out
    a test split.
  * a HARD-NEGATIVE set of related-but-distinct protein families (e.g. other
    GNAT acetyltransferases for ectA, other class-III aminotransferases for
    ectB) — used to calibrate HMM score cutoffs and as specificity controls
    in the benchmark.

Usage:  python 01_fetch_refs.py [--families families.yaml] [--out refs]
Output: refs/<family>.positive.faa
        refs/<family>.negative.faa
        refs/manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import time
from datetime import date
from pathlib import Path

import requests
import yaml

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_INFO = "https://rest.uniprot.org/utils/release"
BATCH_SIZE = 500
SLEEP_BETWEEN_PAGES = 0.5


def load_families(path: Path) -> list[dict]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return data["families"]


def get_uniprot_release() -> str:
    try:
        r = requests.get(UNIPROT_INFO, timeout=10)
        r.raise_for_status()
        return r.json().get("releaseNumber", "unknown")
    except Exception:
        return "unknown"


def fetch_all_sequences(query: str, max_seqs: int | None = None) -> str:
    """
    Paginate through UniProt results for a query and return combined FASTA
    text. If max_seqs is given, stop once at least that many sequences have
    been fetched (some negative-family queries match hundreds of thousands
    of UniProt entries; a representative sample is enough for cutoff
    calibration / benchmarking).
    """
    params = {"query": query, "format": "fasta", "size": BATCH_SIZE}
    all_fasta = []
    url = UNIPROT_API
    page = 1
    n_fetched = 0

    while url:
        r = requests.get(url, params=params if page == 1 else None, timeout=60)
        r.raise_for_status()
        all_fasta.append(r.text)
        n_fetched += r.text.count(">")

        if max_seqs is not None and n_fetched >= max_seqs:
            break

        link_header = r.headers.get("Link", "")
        next_url = None
        for part in link_header.split(","):
            part = part.strip()
            if 'rel="next"' in part:
                next_url = part.split(";")[0].strip().strip("<>")
                break
        url = next_url
        page += 1
        if url:
            time.sleep(SLEEP_BETWEEN_PAGES)

    return "".join(all_fasta)


def retag_fasta(fasta_text: str, tag: str) -> str:
    """Retag FASTA headers to: >tag|UniProtID|OrganismName"""
    lines = fasta_text.strip().split("\n")
    retagged = []
    for line in lines:
        if not line:
            continue
        if line.startswith(">"):
            parts = line[1:].split("|")
            uniprot_id = parts[1] if len(parts) >= 2 else parts[0].strip()
            os_start = line.find("OS=")
            ox_start = line.find("OX=")
            if os_start != -1 and ox_start != -1:
                organism = line[os_start + 3:ox_start].strip().replace(" ", "_")
            elif os_start != -1:
                organism = line[os_start + 3:].split()[0]
            else:
                organism = "unknown"
            retagged.append(f">{tag}|{uniprot_id}|{organism}")
        else:
            retagged.append(line)
    return "\n".join(retagged) + "\n"


def count_seqs(fasta_text: str) -> int:
    return fasta_text.count("\n>") + (1 if fasta_text.lstrip().startswith(">") else 0)


def fetch_set(query: str, tag: str, label: str, max_seqs: int | None = None) -> tuple[str, int]:
    print(f"  {label} query: {query}")
    try:
        fasta_raw = fetch_all_sequences(query, max_seqs=max_seqs)
        if not fasta_raw.strip():
            print(f"  WARNING: no sequences returned for {label} — check query")
            return "", 0
        fasta_retagged = retag_fasta(fasta_raw, tag)
        n = count_seqs(fasta_retagged)
        print(f"  Retrieved: {n} sequences")
        return fasta_retagged, n
    except requests.RequestException as e:
        print(f"  ERROR fetching {label}: {e}")
        return "", 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--out", type=Path, default=Path("refs"))
    ap.add_argument("--max-positive", type=int, default=None,
                    help="Cap on number of positive sequences fetched per family.")
    ap.add_argument("--max-negative", type=int, default=1000,
                    help="Cap on number of hard-negative sequences fetched per family.")
    args = ap.parse_args()

    families = load_families(args.families)
    args.out.mkdir(parents=True, exist_ok=True)
    fetch_date = date.today().isoformat()

    print("Fetching UniProt release info...")
    uniprot_release = get_uniprot_release()
    print(f"UniProt release: {uniprot_release}")
    print(f"Fetch date:      {fetch_date}\n")

    manifest_rows = []

    for fam in families:
        name = fam["name"]
        pos_query = fam["positive_query"]
        neg_query = fam["negative_query"]
        description = fam.get("description", "").strip()

        print(f"[{name}]")

        pos_fasta, n_pos = fetch_set(pos_query, name, "positive", max_seqs=args.max_positive)
        (args.out / f"{name}.positive.faa").write_text(pos_fasta)

        neg_fasta, n_neg = fetch_set(neg_query, f"{name}_neg", "negative", max_seqs=args.max_negative)
        (args.out / f"{name}.negative.faa").write_text(neg_fasta)

        manifest_rows.append({
            "family": name,
            "n_positive": n_pos,
            "n_negative": n_neg,
            "uniprot_release": uniprot_release,
            "date_fetched": fetch_date,
            "positive_query": pos_query,
            "negative_query": neg_query,
            "description": description,
        })
        print()

    manifest_path = args.out / "manifest.tsv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    total_pos = sum(r["n_positive"] for r in manifest_rows)
    total_neg = sum(r["n_negative"] for r in manifest_rows)
    print(f"Done. {total_pos} positive / {total_neg} negative sequences across "
          f"{len(manifest_rows)} families.")
    print(f"Manifest written to {manifest_path}")


if __name__ == "__main__":
    main()
