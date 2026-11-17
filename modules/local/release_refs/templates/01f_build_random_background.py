#!/usr/bin/env python3
"""
01f_build_random_background.py — shared random bacterial background for the benchmark

Each family's own negatives are hard negatives: proteins sharing its Pfam
domain. They are what cutoffs must be calibrated against, but scored alone
they only say how well a family is told apart from its closest relatives.
They say nothing about the reads that make up nearly all of a real
metagenome: proteins unrelated to anything in the panel. This step builds
that second, shared negative set, used only by the benchmark (never by
calibration, never in the DIAMOND or HMM db):

  1. Sample --n random bacterial proteins from UniProt that carry none of
     the Pfam domains named anywhere in the families file (positive,
     negative and hmm-negative queries, negative_pfam, fusion_marker_pfam,
     pfam_model). Uses 01_fetch_refs.py's sampler, so a population this
     size (~75 million) is drawn by stratified sampling and cached per
     UniProt release.
  2. Drop any sampled protein >= --threshold % identical (DIAMOND) to any
     family's positive set, the same rule 01b_check_negative_purity.py uses
     for hard negatives. This catches panel genes the domain exclusion
     misses: genes queried by symbol whose Pfam domain isn't listed (e.g.
     the marker set), or entries with no Pfam annotation at all.

09_simulate_reads.py simulates reads from the result as family
"random_bg"; 11_compute_metrics.py counts every call on those reads as an
off-target false positive for the family it was called as.

Outputs (in --refs):
  random_background.faa           kept background proteins (">random_bg|ACC|Organism")
  random_background.flagged.faa   dropped as too similar to a panel positive
  random_background_manifest.tsv  query, release, date, counts

Usage:
  python 01f_build_random_background.py --refs refs --families families.yaml [--n 5000]
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import random
import re
import shutil
import subprocess
import tempfile
from datetime import date
from pathlib import Path

import yaml

RANDOM_BG_FAMILY = "random_bg"           # must match 09_simulate_reads.py / 11
RANDOM_BG_FASTA = "random_background.faa"
PFAM_FIELDS = ("positive_query", "negative_query", "hmm_negative_query",
               "negative_pfam", "fusion_marker_pfam", "pfam_model")

HERE = Path(__file__).resolve().parent


def load_fetch_module():
    """01_fetch_refs.py isn't importable by name (leading digit)."""
    spec = importlib.util.spec_from_file_location("fetch_refs", HERE / "01_fetch_refs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def panel_pfam_domains(families: list[dict]) -> list[str]:
    domains: set[str] = set()
    for fam in families:
        for field in PFAM_FIELDS:
            value = fam.get(field)
            if value:
                domains.update(re.findall(r"PF\d{5}", str(value)))
    return sorted(domains)


def background_query(domains: list[str]) -> str:
    if not domains:
        return "taxonomy_id:2"
    if len(domains) > 100:
        # UniProt rejects more than 100 OR conditions in one query.
        raise SystemExit(f"ERROR: {len(domains)} panel Pfam domains, more than UniProt's "
                         f"100-condition OR limit; split the exclusion")
    return "taxonomy_id:2 NOT (" + " OR ".join(f"xref:pfam-{d}" for d in domains) + ")"


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records, header, seq = [], None, []
    for line in text.splitlines():
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


def flag_similar_to_positives(records: list[tuple[str, str]], positive_paths: list[Path],
                              threshold: float, threads: int) -> set[str]:
    """Headers of background proteins >= threshold % identical to any positive."""
    if not shutil.which("diamond"):
        raise SystemExit("ERROR: diamond not on PATH (needed for the background purity check)")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        pos_all = tmp_dir / "positives.faa"
        with open(pos_all, "w") as out:
            for path in positive_paths:
                out.write(path.read_text())
        bg = tmp_dir / "background.faa"
        write_fasta(bg, records)
        subprocess.run(["diamond", "makedb", "--in", str(pos_all), "--db",
                        str(tmp_dir / "pos.dmnd"), "--quiet"], check=True)
        hits = tmp_dir / "hits.tsv"
        subprocess.run(["diamond", "blastp", "--query", str(bg), "--db", str(tmp_dir / "pos.dmnd"),
                        "--out", str(hits), "--outfmt", "6", "qseqid", "pident",
                        "--max-target-seqs", "1", "--threads", str(threads), "--quiet"],
                       check=True)
        flagged = set()
        for line in hits.read_text().splitlines():
            qseqid, pident = line.split("\t")
            if float(pident) >= threshold:
                flagged.add(qseqid)
    # DIAMOND reports qseqid up to the first whitespace; headers have none.
    return flagged


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--n", type=int, default=5000,
                    help="Background proteins to sample before the purity check")
    ap.add_argument("--threshold", type=float, default=70.0,
                    help="Drop background proteins at or above this %% identity "
                         "to any panel positive (01b's rule)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--cache-dir", type=Path, default=None,
                    help="Sampling cache (default: 01_fetch_refs.py's)")
    args = ap.parse_args()

    fetch = load_fetch_module()
    with open(args.families) as fh:
        families = yaml.safe_load(fh)["families"]

    fetch.UNIPROT_RELEASE = fetch.get_uniprot_release()
    fetch.CACHE_DIR = args.cache_dir or fetch.DEFAULT_CACHE_DIR
    domains = panel_pfam_domains(families)
    query = background_query(domains)
    print(f"UniProt release: {fetch.UNIPROT_RELEASE}")
    print(f"Excluding {len(domains)} panel Pfam domains")

    fasta, n_sampled, _ = fetch.fetch_set(
        query, RANDOM_BG_FAMILY, "random background", max_seqs=args.n,
        rng=random.Random(f"{fetch.FETCH_RANDOM_SEED}:{RANDOM_BG_FAMILY}"))
    records = parse_fasta(fasta)

    positive_paths = [args.refs / f"{fam['name']}.positive.faa" for fam in families]
    missing = [p.name for p in positive_paths if not p.exists()]
    if missing:
        raise SystemExit(f"ERROR: positives missing ({', '.join(missing[:5])}...): run steps 01-01e first")
    flagged_ids = flag_similar_to_positives(records, positive_paths, args.threshold, args.threads)
    kept = [r for r in records if r[0] not in flagged_ids]
    flagged = [r for r in records if r[0] in flagged_ids]

    write_fasta(args.refs / RANDOM_BG_FASTA, kept)
    write_fasta(args.refs / "random_background.flagged.faa", flagged)
    with open(args.refs / "random_background_manifest.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["query", "uniprot_release", "date_fetched", "n_sampled",
                         "n_flagged_similar_to_positive", "n_kept", "identity_threshold"])
        writer.writerow([query, fetch.UNIPROT_RELEASE, date.today().isoformat(), n_sampled,
                         len(flagged), len(kept), args.threshold])
    print(f"Random background: {n_sampled} sampled, {len(flagged)} dropped as >= "
          f"{args.threshold:.0f}% identical to a panel positive, {len(kept)} kept -> "
          f"{args.refs / RANDOM_BG_FASTA}")


if __name__ == "__main__":
    main()
