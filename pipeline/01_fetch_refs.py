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
import random
import re
import time
from datetime import date
from pathlib import Path

import requests
import yaml

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_INFO = "https://rest.uniprot.org/utils/release"
BATCH_SIZE = 500
ACCESSION_BATCH_SIZE = 100   # UniProt REST hard-caps OR queries at 100 conditions (confirmed: "Too many OR
                              # conditions in query. Maximum allowed is 100." for 200) -- not a URL-length limit
SLEEP_BETWEEN_PAGES = 0.5
FETCH_RANDOM_SEED = 42


def next_page_url(link_header: str) -> str | None:
    """Extract the rel="next" URL from a UniProt Link header.

    NOT a naive split(",") -- a request with multiple comma-separated
    `fields` (e.g. "accession,xref_pfam", used by fetch_pfam_domains below)
    produces a Link header whose URL itself contains a comma, which a plain
    split(",") chops mid-URL. Match each "<url>; rel=..." segment instead,
    since the only commas that matter are the ones separating whole link
    entries, immediately after a closing '>'."""
    for url, rel in re.findall(r'<([^>]+)>\s*;\s*rel="([^"]+)"', link_header):
        if rel == "next":
            return url
    return None


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


def fetch_all_sequences(query: str) -> str:
    """Paginate through ALL UniProt results for a query and return combined
    FASTA text. No cap -- see resolve_accessions() for how a family-level
    cap (max_negative / max_positive_override) is applied instead, via
    accession-level random sampling rather than early-stopping pagination."""
    params = {"query": query, "format": "fasta", "size": BATCH_SIZE}
    all_fasta = []
    url = UNIPROT_API
    page = 1

    while url:
        r = requests.get(url, params=params if page == 1 else None, timeout=60)
        r.raise_for_status()
        all_fasta.append(r.text)

        url = next_page_url(r.headers.get("Link", ""))
        page += 1
        if url:
            time.sleep(SLEEP_BETWEEN_PAGES)

    return "".join(all_fasta)


def fetch_all_accessions(query: str, sort: str | None = None, stop_after: int | None = None) -> list[str]:
    """Paginate through matching accessions (lightweight: accession field
    only, no sequence data). With stop_after=None, fetches the ENTIRE
    population -- used as the shared basis for both the FASTA and (for
    fusion_partner families) Pfam-domain fetches, so they agree on exactly
    the same sequences, for queries small enough that's practical. With
    stop_after set, stops once at least that many accessions are
    collected -- see resolve_accessions() for why/when that matters."""
    params = {"query": query, "format": "tsv", "fields": "accession", "size": BATCH_SIZE}
    if sort:
        params["sort"] = sort
    accessions: list[str] = []
    url = UNIPROT_API
    page = 1

    while url:
        r = requests.get(url, params=params if page == 1 else None, timeout=60)
        r.raise_for_status()
        accessions.extend(line for line in r.text.strip().split("\n")[1:] if line)

        if stop_after is not None and len(accessions) >= stop_after:
            break

        url = next_page_url(r.headers.get("Link", ""))
        page += 1
        if url:
            time.sleep(SLEEP_BETWEEN_PAGES)

    return accessions


def fetch_fasta_by_accessions(accessions: list[str]) -> str:
    """Fetch FASTA for a specific, already-chosen list of accessions, in
    ACCESSION_BATCH_SIZE-sized `accession:(A OR B OR ...)` queries (well
    under typical URL length limits)."""
    chunks = []
    for i in range(0, len(accessions), ACCESSION_BATCH_SIZE):
        batch = accessions[i:i + ACCESSION_BATCH_SIZE]
        query = "accession:(" + " OR ".join(batch) + ")"
        r = requests.get(UNIPROT_API, params={"query": query, "format": "fasta", "size": ACCESSION_BATCH_SIZE},
                          timeout=60)
        r.raise_for_status()
        chunks.append(r.text)
        if i + ACCESSION_BATCH_SIZE < len(accessions):
            time.sleep(SLEEP_BETWEEN_PAGES)
    return "".join(chunks)


OVERSAMPLE_FACTOR = 5   # how many candidates to page through per accession actually wanted, before sampling


def resolve_accessions(query: str, max_seqs: int | None, rng: random.Random) -> list[str]:
    """All matching accessions for query, or -- if max_seqs is set -- a
    random sample of that size drawn from a bounded oversample.

    UniProt's default result order (no explicit sort param) is NOT
    representative: it's heavily biased toward reviewed/well-studied model
    organisms -- confirmed directly (a plain cspA-style Pfam query put E.
    coli/B. subtilis in 8 of its first 10 hits). Simply stopping pagination
    early once max_seqs results have been SEEN (the original approach here)
    silently inherited that bias, which would undercut exactly the
    "diverse gut taxa" goal a capped family like cspA exists for -- and,
    since every family's negative fetch is capped by --max-negative
    (default 1000) too, was quietly affecting hard-negative diversity
    panel-wide, not just cspA's positives.

    Enumerating the ENTIRE population first (the first fix attempted here)
    is unbiased but was found to be impractically slow for a broad Pfam
    family's negative pool: ectA's alone (PF00583) has 311,116 UniProt
    members, which at 500 accessions/page with rate-limit-friendly
    pagination is 15-25+ minutes just to LIST -- before any sequences are
    even fetched -- multiplied across every broad-domain family in the
    panel. Confirmed via direct testing that `sort=accession asc` alone
    (no exhaustive enumeration needed) already fixes the actual bias
    cheaply: a 200-accession sample under that sort spanned dozens of
    genera with no single genus over 12/200, vastly better than the
    default order's 8/10 in two species, without needing to know the
    population size at all. So: page through accession-sorted results
    until OVERSAMPLE_FACTOR * max_seqs candidates are collected (bounded,
    cheap, independent of true population size), then randomly sample
    max_seqs from that oversample -- unbiased order + genuine random draw,
    without ever needing to enumerate a 300k-member pool in full."""
    if max_seqs is None:
        return fetch_all_accessions(query)
    accessions = fetch_all_accessions(query, sort="accession asc", stop_after=max_seqs * OVERSAMPLE_FACTOR)
    if len(accessions) > max_seqs:
        accessions = rng.sample(accessions, max_seqs)
    return accessions


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


def diversity_stats(fasta_text: str) -> tuple[int, int]:
    """(n_distinct_organisms, n_distinct_genera) from a retagged FASTA's own
    '>tag|accession|Organism_name' headers -- free (no extra network calls,
    the organism is already embedded by retag_fasta) and exactly the kind
    of check that would have caught the default-sort-order bias fixed
    above at build time instead of by manual investigation after the fact.
    Genus = organism's first underscore-separated token, same convention
    03_split_train_test.py's genus_of() uses for taxonomy-mode splitting."""
    organisms, genera = set(), set()
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            organism = parts[2] if len(parts) >= 3 else "unknown"
            organisms.add(organism)
            genera.add(organism.split("_")[0])
    return len(organisms), len(genera)


def count_seqs(fasta_text: str) -> int:
    return fasta_text.count("\n>") + (1 if fasta_text.lstrip().startswith(">") else 0)


def extract_accessions_from_fasta(fasta_text: str) -> list[str]:
    """Pull the UniProt accession out of each raw '>sp|ACCESSION|...' /
    '>tr|ACCESSION|...' header -- same parsing rule retag_fasta uses."""
    accessions = []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            accessions.append(parts[1] if len(parts) >= 2 else parts[0].strip())
    return accessions


def fetch_set(query: str, tag: str, label: str, max_seqs: int | None = None,
              rng: random.Random | None = None) -> tuple[str, int, list[str]]:
    """Returns (retagged FASTA text, n sequences, accessions actually used) --
    the accession list is returned so a fusion_partner family's domain-
    evidence fetch (see fetch_pfam_domains) can reuse the exact same set,
    whether or not this fetch happened to be capped/sampled."""
    print(f"  {label} query: {query}")
    try:
        if max_seqs is not None:
            accessions = resolve_accessions(query, max_seqs, rng or random.Random(FETCH_RANDOM_SEED))
            if not accessions:
                print(f"  WARNING: no sequences returned for {label} — check query")
                return "", 0, []
            fasta_raw = fetch_fasta_by_accessions(accessions)
        else:
            fasta_raw = fetch_all_sequences(query)
            if not fasta_raw.strip():
                print(f"  WARNING: no sequences returned for {label} — check query")
                return "", 0, []
            accessions = extract_accessions_from_fasta(fasta_raw)
        fasta_retagged = retag_fasta(fasta_raw, tag)
        n = count_seqs(fasta_retagged)
        print(f"  Retrieved: {n} sequences")
        return fasta_retagged, n, accessions
    except requests.RequestException as e:
        print(f"  ERROR fetching {label}: {e}")
        return "", 0, []


def fetch_pfam_domains(accessions: list[str]) -> dict[str, list[str]]:
    """
    Fetch each accession's Pfam domain-accession list (UniProt's own
    curated annotation, not something we compute), for an EXPLICIT list of
    accessions -- the same ones actually fetched into positive.faa by
    fetch_set above, so domain evidence never disagrees with what's really
    in the reference set (matters once a family is random-sampled via
    max_positive_override). Used only for families that declare
    `fusion_partner` in families.yaml, to give
    01c_check_length_outliers.py real domain evidence (not just length) for
    telling a genuine fused ORF apart from an unrelated length outlier --
    see that script's docstring.
    """
    domains: dict[str, list[str]] = {}
    for i in range(0, len(accessions), ACCESSION_BATCH_SIZE):
        batch = accessions[i:i + ACCESSION_BATCH_SIZE]
        query = "accession:(" + " OR ".join(batch) + ")"
        r = requests.get(UNIPROT_API,
                          params={"query": query, "format": "tsv", "fields": "accession,xref_pfam",
                                  "size": ACCESSION_BATCH_SIZE},
                          timeout=60)
        r.raise_for_status()
        for line in r.text.strip().split("\n")[1:]:  # skip header
            if not line:
                continue
            parts = line.split("\t")
            accession = parts[0]
            pfam = parts[1].strip(";").split(";") if len(parts) > 1 and parts[1] else []
            domains[accession] = pfam
        if i + ACCESSION_BATCH_SIZE < len(accessions):
            time.sleep(SLEEP_BETWEEN_PAGES)

    return domains


def check_no_duplicate_pfam_models(families: list[dict]) -> None:
    """
    Fail fast if two families share a `pfam_model` accession.

    05_build_hmms.sh fetches the literal Pfam-A HMM for a pfam_model family
    and only rewrites its NAME/ACC lines -- the actual model weights are
    untouched, so two families adopting the SAME accession end up with
    byte-identical HMMs. hmmscan then can't discriminate between them at
    all: every read scores exactly the same against both, and whichever
    sorts first in 07_press_hmms.sh's alphabetical `cat hmms/*.hmm`
    ordering wins 100% of the shared signal via 11_compute_metrics.py's
    best-hit tie-break, leaving the other(s) at exactly zero HMM recall.
    Confirmed in production (v6 benchmark): trkH/ktrB/ktrD all adopted
    pfam_model: PF02386 -- ktrB absorbed every tied hit, trkH and ktrD
    both scored 0 true positives despite real, non-trivial DIAMOND recall
    on the same reads. Fixed by removing pfam_model from all three (see
    families.yaml). If a future family group is tempted to make the same
    choice: at most one member of a shared-Pfam-accession sibling group
    (see families.yaml's proX/opuAC/opuBC/opuCC precedent) may use
    pfam_model for that accession; the rest must use a custom-built model.
    """
    by_accession: dict[str, list[str]] = {}
    for fam in families:
        accession = fam.get("pfam_model")
        if accession:
            by_accession.setdefault(accession, []).append(fam["name"])
    conflicts = {acc: names for acc, names in by_accession.items() if len(names) > 1}
    if conflicts:
        lines = [f"  {acc}: {', '.join(names)}" for acc, names in conflicts.items()]
        raise SystemExit(
            "ERROR: multiple families.yaml entries share a pfam_model accession -- "
            "this makes their HMMs byte-identical and breaks HMM-based detection "
            "for all but one of them (see check_no_duplicate_pfam_models docstring "
            "for the confirmed failure mode):\n" + "\n".join(lines)
        )


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
    check_no_duplicate_pfam_models(families)
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
        # Per-family cap, for a family whose positive_query is deliberately
        # anchored on a broad Pfam accession rather than a gene symbol (e.g.
        # cspA-family: xref:pfam-PF00313 alone matches ~80k UniProt proteins)
        # -- overrides --max-positive for that family only.
        max_pos = fam.get("max_positive_override", args.max_positive)

        print(f"[{name}]")

        # Per-family, per-label RNG so a random-sampled fetch is
        # reproducible across reruns but independent of every other
        # family/label's own draws (same reasoning as
        # 03_split_train_test.py's per-family seeding).
        pos_fasta, n_pos, pos_accessions = fetch_set(
            pos_query, name, "positive", max_seqs=max_pos,
            rng=random.Random(f"{FETCH_RANDOM_SEED}:{name}:positive"))
        (args.out / f"{name}.positive.faa").write_text(pos_fasta)
        n_pos_orgs, n_pos_genera = diversity_stats(pos_fasta)
        print(f"  Diversity: {n_pos_orgs} organisms / {n_pos_genera} genera")

        neg_fasta, n_neg, _ = fetch_set(
            neg_query, f"{name}_neg", "negative", max_seqs=args.max_negative,
            rng=random.Random(f"{FETCH_RANDOM_SEED}:{name}:negative"))
        (args.out / f"{name}.negative.faa").write_text(neg_fasta)
        n_neg_orgs, n_neg_genera = diversity_stats(neg_fasta)
        print(f"  Diversity: {n_neg_orgs} organisms / {n_neg_genera} genera")

        # Fusion-partner families (e.g. mrpA/mrpB, see families.yaml) need
        # real Pfam domain evidence, not just length, to tell a genuine
        # fused ORF apart from an unrelated length outlier in step 01c --
        # fetch it now while we have network access to UniProt, for the
        # exact same accessions that ended up in positive.faa above.
        if fam.get("fusion_partner"):
            domains = fetch_pfam_domains(pos_accessions)
            domains_path = args.out / f"{name}.positive.domains.tsv"
            with open(domains_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(["accession", "pfam_domains"])
                for accession, pfam in domains.items():
                    writer.writerow([accession, ";".join(pfam)])
            print(f"  Pfam domain evidence for {len(domains)} accessions -> {domains_path}")

        manifest_rows.append({
            "family": name,
            "n_positive": n_pos,
            "n_positive_organisms": n_pos_orgs,
            "n_positive_genera": n_pos_genera,
            "n_negative": n_neg,
            "n_negative_organisms": n_neg_orgs,
            "n_negative_genera": n_neg_genera,
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
