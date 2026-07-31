#!/usr/bin/env python3
"""
01e_add_refseq_positives.py -- merge additional real positives from NCBI
RefSeq into a family's positive pool, for families whose UniProt gene-symbol
population is already fully exhausted (fetched count == UniProt's own total,
no fetch cap left to raise) but where RefSeq's independently-annotated
bacterial protein set turned out to add real, non-redundant, correctly-
annotated sequences.

Only families that declare `refseq_gene_symbols: [...]` in families.yaml are
fetched here -- this is NOT run unconditionally for every family. A bare gene
symbol can collide with a completely unrelated gene in a different organism:
confirmed during evaluation that "MrpC" is also a Proteus fimbrial usher gene
and a Myxococcus CRP/FNR transcription factor, and "PhaF" (mrpF's own
documented alias) is also a common polyhydroxyalkanoate-granule protein name.
Each enabled family was checked individually first (CD-HIT-2D novelty check
against the existing pool + a manual product-description spot-check) -- see
docs/CHANGELOG.md and each family's own families.yaml comment.

Sampling: NCBI's esearch default result order is not randomized (the same
kind of order bias 01_fetch_refs.py already had to fix for UniProt -- see
issue #4 in docs/CHANGELOG.md), so this fetches the FULL matching UID list
first (lightweight -- UIDs only, no sequence data) and draws a reproducible
random sample from it, the same oversample-free approach 01_fetch_refs.py
uses once the full population is already in hand.

Entries whose RefSeq description is explicitly marked ", partial" are
skipped at fetch time (a free, unambiguous fragment signal) -- everything
else still goes through the normal length-outlier filter downstream.

Usage:
  python 01e_add_refseq_positives.py --refs refs --families families.yaml
Output:
  refs/<family>.positive.faa           (RefSeq sequences appended in place)
  refs/refseq_positives_manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import random
import re
import time
from datetime import date
from pathlib import Path

import requests
import yaml

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_BATCH_SIZE = 100000
EFETCH_BATCH_SIZE = 200
SLEEP_BETWEEN_REQUESTS = 0.4
REFSEQ_RANDOM_SEED = 42


def load_families(path: Path) -> list[dict]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return data["families"]


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


def append_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "a") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def existing_ids(records: list[tuple[str, str]]) -> set[str]:
    ids = set()
    for header, _ in records:
        parts = header.split("|")
        ids.add(parts[1] if len(parts) >= 2 else header)
    return ids


def build_query(gene_symbols: list[str]) -> str:
    gene_clause = " OR ".join(f"{sym}[Gene Name]" for sym in gene_symbols)
    return f"({gene_clause}) AND bacteria[Organism] AND srcdb_refseq[Properties]"


MAX_RETRIES = 5


def _request_with_retry(method: str, url: str, **kwargs) -> requests.Response:
    """NCBI's eutils endpoints are known to return transient 5xx/429s under
    normal load (observed directly: a 502 mid-run here), independent of
    anything this script does wrong -- retry with exponential backoff rather
    than letting one flaky response abort an otherwise-successful family."""
    last_exc = None
    for attempt in range(MAX_RETRIES):
        try:
            r = requests.request(method, url, **kwargs)
            r.raise_for_status()
            return r
        except requests.exceptions.RequestException as exc:
            last_exc = exc
            if attempt < MAX_RETRIES - 1:
                wait = SLEEP_BETWEEN_REQUESTS * (2 ** attempt)
                print(f"    (transient error: {exc} -- retrying in {wait:.1f}s)")
                time.sleep(wait)
    raise last_exc


def esearch_all_uids(query: str) -> list[str]:
    """All matching UIDs for query -- lightweight (UIDs only), so fetching
    the full list (rather than a capped/paginated subset) is cheap even for
    a broad query, and lets sampling below draw an unbiased random subset
    instead of trusting esearch's own (non-random) default ordering."""
    r = _request_with_retry("GET", f"{EUTILS}/esearch.fcgi",
                             params={"db": "protein", "term": query, "retmax": ESEARCH_BATCH_SIZE},
                             timeout=60)
    return re.findall(r"<Id>(\d+)</Id>", r.text)


def efetch_fasta(uids: list[str]) -> str:
    chunks = []
    for i in range(0, len(uids), EFETCH_BATCH_SIZE):
        batch = uids[i:i + EFETCH_BATCH_SIZE]
        r = _request_with_retry("POST", f"{EUTILS}/efetch.fcgi",
                                 data={"db": "protein", "id": ",".join(batch),
                                       "rettype": "fasta", "retmode": "text"},
                                 timeout=60)
        chunks.append(r.text)
        if i + EFETCH_BATCH_SIZE < len(uids):
            time.sleep(SLEEP_BETWEEN_REQUESTS)
    return "".join(chunks)


HEADER_RE = re.compile(r"^(\S+) (.+?) \[([^\]]+)\]\s*$")


def parse_refseq_fasta(fasta_text: str) -> list[tuple[str, str, str]]:
    """Returns (accession, organism, sequence) for each non-partial record."""
    records = []
    accession = organism = None
    partial = False
    chunks: list[str] = []

    def flush():
        if accession is not None and not partial and chunks:
            records.append((accession, organism, "".join(chunks)))

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            flush()
            chunks = []
            m = HEADER_RE.match(line[1:].strip())
            if m:
                accession, product, organism = m.groups()
                partial = product.rstrip().endswith(", partial")
                organism = organism.replace(" ", "_")
            else:
                accession, organism, partial = line[1:].split()[0], "unknown", False
        elif line.strip():
            chunks.append(line.strip())
    flush()
    return records


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--max-per-family", type=int, default=1500,
                     help="Max RefSeq sequences to sample per family (before dedup)")
    args = ap.parse_args()

    fam_defs = [f for f in load_families(args.families) if f.get("refseq_gene_symbols")]
    if not fam_defs:
        print("No family declares refseq_gene_symbols -- nothing to fetch.")
        return

    # NCBI has no discrete, citable "release number" for RefSeq the way
    # UniProt does (it's continuously updated) -- a fetch date is the
    # practical provenance signal here, same role date_fetched already
    # plays in 01_fetch_refs.py's own manifest.
    date_fetched = date.today().isoformat()
    print(f"Fetch date: {date_fetched}")

    manifest_rows = []
    for fam in fam_defs:
        name = fam["name"]
        pos_path = args.refs / f"{name}.positive.faa"
        if not pos_path.exists():
            print(f"[{name}] SKIP: {pos_path} does not exist yet -- run 01_fetch_refs.py first")
            continue

        query = build_query(fam["refseq_gene_symbols"])
        print(f"[{name}] RefSeq query: {query}")
        uids = esearch_all_uids(query)
        n_total = len(uids)

        rng = random.Random(f"{REFSEQ_RANDOM_SEED}:{name}:refseq")
        sample_uids = uids if n_total <= args.max_per_family else rng.sample(uids, args.max_per_family)
        print(f"[{name}] RefSeq population: {n_total}, sampling {len(sample_uids)}")

        fasta_text = efetch_fasta(sample_uids)
        refseq_records = parse_refseq_fasta(fasta_text)
        n_partial_skipped = len(sample_uids) - len(refseq_records)

        current_ids = existing_ids(parse_fasta(pos_path))
        to_add = []
        n_dupe = 0
        for accession, organism, seq in refseq_records:
            if accession in current_ids:
                n_dupe += 1
                continue
            to_add.append((f"{name}_refseq|{accession}|{organism}", seq))
            current_ids.add(accession)

        append_fasta(pos_path, to_add)
        print(f"[{name}] {len(refseq_records)} fetched (non-partial): "
              f"{len(to_add)} appended, {n_dupe} already present, "
              f"{n_partial_skipped} partial skipped")

        manifest_rows.append({
            "family": name,
            "refseq_population": n_total,
            "n_sampled": len(sample_uids),
            "n_partial_skipped": n_partial_skipped,
            "n_appended": len(to_add),
            "n_already_present": n_dupe,
            "date_fetched": date_fetched,
        })
        time.sleep(SLEEP_BETWEEN_REQUESTS)

    if not manifest_rows:
        print("No families with refseq_gene_symbols had an existing positive.faa to merge into.")
        return

    manifest_path = args.refs / "refseq_positives_manifest.tsv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    n_total_added = sum(r["n_appended"] for r in manifest_rows)
    print(f"\nDone. {n_total_added} RefSeq positive sequences merged across {len(manifest_rows)} families.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
