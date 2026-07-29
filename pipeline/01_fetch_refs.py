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
        refs/<family>.positive.domains.tsv, refs/<family>.negative.domains.tsv
                                            (fusion_partner families only)
        refs/<famA>_<famB>.dedicated_fusion_fetch.faa
                                            (declared fusion_partner pairs
                                            only -- see fetch_fusion_pair_
                                            candidates)
        refs/manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import statistics
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


def fusion_pairs(families: list[dict]) -> list[tuple[str, str]]:
    """Unique, sorted fusion_partner pairs declared in families.yaml -- same
    dedup logic as 08d_build_fusion_refs.py's load_fusion_pairs (kept
    independent rather than shared, since the two scripts have no import
    relationship)."""
    fam_by_name = {fam["name"]: fam for fam in families}
    seen = set()
    pairs = []
    for fam in families:
        partner = fam.get("fusion_partner")
        if not partner or partner not in fam_by_name:
            continue
        key = frozenset((fam["name"], partner))
        if key in seen:
            continue
        seen.add(key)
        pairs.append(tuple(sorted((fam["name"], partner))))
    return pairs


def load_positive_domain_fractions(out_dir: Path, family: str) -> dict[str, float]:
    """pfam_accession -> fraction of `family`'s own gene-symbol-matched
    positive fetch that carries it, from the domains.tsv this same run just
    wrote (fetch_pfam_domains on pos_accessions, see main()). Used to tell
    a domain that's genuinely PART OF a family's own standalone
    architecture (fraction near 1.0) apart from one that's foreign to it
    (fraction near 0) -- see fetch_fusion_pair_candidates.

    Restricted to the "typical-length" subset of the raw positive fetch
    (within [median/1.5, median*1.5], the same window 01c_check_length_
    outliers.py itself uses) before computing fractions -- NOT the whole
    raw fetch. A fusion_partner family's own gene-symbol query can already
    include a real, non-trivial fraction of already-tagged fused ORFs
    (confirmed: mrpA's raw positive fetch is ~35% oversized, already-tagged
    fusion entries, e.g. "gene:mrpA" on a 957aa Paenibacillus apis ORF --
    see families.yaml's mrpA entry). Computing fractions over the raw fetch
    directly is circular for exactly this reason: it would count the
    partner's own domain as "not foreign" simply because the raw fetch is
    already partly contaminated with fusions, defeating the point of this
    check. Restricting to typical-length entries first removes that
    circularity."""
    domains_path = out_dir / f"{family}.positive.domains.tsv"
    fasta_path = out_dir / f"{family}.positive.faa"
    if not domains_path.exists() or not fasta_path.exists():
        return {}

    lengths = dict(_parse_fasta_lengths(fasta_path))
    if not lengths:
        return {}
    median_len = statistics.median(lengths.values())
    lower, upper = median_len / 1.5, median_len * 1.5
    typical_accessions = {acc for acc, length in lengths.items() if lower <= length <= upper}

    counts: dict[str, int] = {}
    total = 0
    with open(domains_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["accession"] not in typical_accessions:
                continue
            total += 1
            for pfam in (row["pfam_domains"].split(";") if row["pfam_domains"] else []):
                counts[pfam] = counts.get(pfam, 0) + 1
    return {pfam: n / total for pfam, n in counts.items()} if total else {}


def _parse_fasta_lengths(path: Path) -> list[tuple[str, int]]:
    """[(accession, sequence_length), ...] from a retagged '>tag|accession|
    organism' FASTA -- accession only (not the full header), to match
    domains.tsv's own accession-keyed rows."""
    records = []
    accession, length = None, 0
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if accession is not None:
                records.append((accession, length))
            parts = line[1:].split("|")
            accession = parts[1] if len(parts) >= 2 else parts[0].strip()
            length = 0
        elif line.strip():
            length += len(line.strip())
    if accession is not None:
        records.append((accession, length))
    return records


def fetch_fusion_pair_candidates(families: list[dict], out_dir: Path) -> None:
    """For each declared fusion_partner pair, fetch a comprehensive
    (uncapped) population of genuine fused-ORF candidates directly, via a
    combined-domain query -- rather than relying only on incidentally
    spotting them in each family's own (capped, gene-symbol-exclusion-based)
    positive/negative fetch, which is what mrpA/mrpB/otsA/otsB relied on
    before this existed. See 01c_check_length_outliers.py's module
    docstring for why incidental discovery alone isn't reliable (it's still
    kept as a safety net -- this is additive, not a replacement).

    Query construction: for pair (A, B), we need one domain that's part of
    A's own standalone architecture but NOT B's, and one that's part of
    B's own standalone architecture but NOT A's -- ANDing those two
    isolates genuine fusions without also matching an ordinary standalone
    member of either family. Candidate domains are each family's own
    fusion_marker_pfam and negative_pfam; which candidate actually has this
    "self-intrinsic to A, foreign to B" property is NOT something families
    .yaml declares directly and can't be assumed from field names alone --
    confirmed the hard way: mrpA and mrpB both declare fusion_marker_pfam
    PF13244, and a naive substitution picked mrpA's own negative_pfam
    (PF00361) as the second domain, but PF00361 is ALSO part of mrpA's own
    standalone architecture (same as PF13244), so "PF13244 AND PF00361"
    just re-matched ~8k ordinary standalone mrpA orthologs instead of
    fusions specifically.

    So this is verified empirically instead, using data already fetched
    this same run: each family's own positive.domains.tsv (gene-symbol-
    matched, i.e. confirmed genuine members), restricted to that family's
    own typical-length subset (see load_positive_domain_fractions) tells us
    what fraction of A's own ordinary members carry a given domain.
    Self-intrinsic to A = that fraction is high; foreign to B = the same
    domain's fraction in B's own typical-length positives is low. Only a
    domain pair passing both checks on both sides is used.

    This still can't be resolved for every declared pair, and that's
    expected, not a bug to chase: mrpA/mrpB fails this check even after
    restricting to typical-length positives, because mrpA is ITSELF
    already a large, multi-domain protein (own median ~800aa) whose own
    size range overlaps substantially with the fused-with-mrpB range
    (~940aa) -- unlike otsA/otsB, where standalone (267aa) and fused
    (900aa+) are night-and-day different. There's no length-based way to
    cleanly separate "ordinary mrpA" from "mrpA already fused" for
    computing this check in the first place, so no domain pair can be
    verified safe. In that case the pair is skipped for the dedicated
    fetch (a printed note, not an error) -- incidental discovery via each
    family's own positive/negative fetch, plus 01c's post-hoc domain-
    evidence pre-filter, remains the (already-validated) mechanism for
    that pair, same as before this dedicated fetch existed.
    """
    SELF_INTRINSIC_MIN = 0.5
    FOREIGN_MAX = 0.1

    fam_by_name = {fam["name"]: fam for fam in families}

    for fam_a, fam_b in fusion_pairs(families):
        a, b = fam_by_name[fam_a], fam_by_name[fam_b]
        out_path = out_dir / f"{fam_a}_{fam_b}.dedicated_fusion_fetch.faa"

        def skip(reason: str) -> None:
            # A previous run may have written this file (e.g. before a
            # families.yaml edit, or before a bug fix here) -- removing it
            # on skip stops 08d_build_fusion_refs.py from silently trusting
            # stale data that this run couldn't actually verify.
            if out_path.exists():
                out_path.unlink()
                reason += f" (removed stale {out_path.name} from a previous run)"
            print(f"[{fam_a}/{fam_b}] SKIP dedicated fusion fetch: {reason}")

        frac_a = load_positive_domain_fractions(out_dir, fam_a)
        frac_b = load_positive_domain_fractions(out_dir, fam_b)
        if not frac_a or not frac_b:
            skip(f"no positive domain evidence for one or both families "
                 f"(need {fam_a}.positive.domains.tsv / {fam_b}.positive.domains.tsv, "
                 f"only written when fusion_partner is declared)")
            continue

        candidates_a = {d for d in (a.get("fusion_marker_pfam"), a.get("negative_pfam")) if d}
        candidates_b = {d for d in (b.get("fusion_marker_pfam"), b.get("negative_pfam")) if d}

        domain_for_a = next((d for d in candidates_a
                              if frac_a.get(d, 0) >= SELF_INTRINSIC_MIN and frac_b.get(d, 0) <= FOREIGN_MAX), None)
        domain_for_b = next((d for d in candidates_b
                              if d != domain_for_a
                              and frac_b.get(d, 0) >= SELF_INTRINSIC_MIN and frac_a.get(d, 0) <= FOREIGN_MAX), None)

        if not domain_for_a or not domain_for_b:
            skip(f"no verified domain pair found where each side is self-intrinsic to one "
                 f"family and foreign to the other (candidates checked: "
                 f"{fam_a}={sorted(candidates_a)}, {fam_b}={sorted(candidates_b)}) -- relying "
                 f"on incidental discovery via each family's own fetch instead")
            continue

        query = f"xref:pfam-{domain_for_a} AND xref:pfam-{domain_for_b} AND taxonomy_id:2"
        print(f"[{fam_a}/{fam_b}] dedicated fusion-candidate query: {query}")
        fasta_raw = fetch_all_sequences(query)
        if not fasta_raw.strip():
            skip(f"no sequences returned for {domain_for_a}/{domain_for_b}")
            continue
        fasta_retagged = retag_fasta(fasta_raw, f"{fam_a}_{fam_b}_dedicated")
        n = count_seqs(fasta_retagged)
        out_path.write_text(fasta_retagged)
        print(f"  {n} candidates -> {out_path}")


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
        # Fusion-partner families: exclude the partner's marker domain from
        # the negative fetch itself, not just after the fact. A genuine
        # fused ORF (bare locus tag, no gene:X symbol) matches this
        # family's own negative_query just as easily as a true hard
        # negative -- confirmed in production for mrpB (~68% of its
        # negative pool) and otsA/otsB (~41-59%), corrupting the
        # median-based length filter in 01c_check_length_outliers.py.
        # Excluding it here means the fetch's own accession-sorted
        # oversample (see resolve_accessions) is no longer wasted on
        # contamination we'd throw away anyway -- the same n=1000 budget
        # now lands entirely on genuine hard negatives. 01c's own
        # domain-evidence pre-filter stays in place as a safety net (e.g.
        # for a fused ORF whose Pfam annotation is stale/lagging at fetch
        # time); see 01c_check_length_outliers.py's module docstring and
        # the dedicated fusion-pair fetch below, which recovers what this
        # exclusion removes as a real detection target rather than just
        # discarding it.
        marker = fam.get("fusion_marker_pfam")
        if fam.get("fusion_partner") and marker:
            neg_query = f"{neg_query} NOT xref:pfam-{marker}"
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

        neg_fasta, n_neg, neg_accessions = fetch_set(
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
        #
        # The NEGATIVE fetch needs the same domain evidence too, not just
        # the positive one: a fused ORF under a bare locus tag (no
        # recognizable gene symbol) matches this family's own negative_query
        # (it carries the Pfam domain the query is anchored on) but not the
        # query's gene-symbol exclusion, so it silently lands in the
        # "negative" pool instead of being caught by the positive-side
        # length+marker check below. Confirmed in production: mrpB's
        # negative pool was ~68% genuine mrpA+mrpB fused-ORF sequences
        # (same PF13244+PF20501+PF04039+PF00361+PF00662 architecture as
        # mrpA's documented fusion cases, under locus tags like
        # "Lokhon_03054") -- see 01c_check_length_outliers.py for what this
        # evidence is used for.
        if fam.get("fusion_partner"):
            domains = fetch_pfam_domains(pos_accessions)
            domains_path = args.out / f"{name}.positive.domains.tsv"
            with open(domains_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(["accession", "pfam_domains"])
                for accession, pfam in domains.items():
                    writer.writerow([accession, ";".join(pfam)])
            print(f"  Pfam domain evidence for {len(domains)} accessions -> {domains_path}")

            neg_domains = fetch_pfam_domains(neg_accessions)
            neg_domains_path = args.out / f"{name}.negative.domains.tsv"
            with open(neg_domains_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(["accession", "pfam_domains"])
                for accession, pfam in neg_domains.items():
                    writer.writerow([accession, ";".join(pfam)])
            print(f"  Pfam domain evidence for {len(neg_domains)} negative accessions -> {neg_domains_path}")

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

    print("Fetching dedicated fusion-pair candidates (comprehensive, uncapped)...")
    fetch_fusion_pair_candidates(families, args.out)
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
