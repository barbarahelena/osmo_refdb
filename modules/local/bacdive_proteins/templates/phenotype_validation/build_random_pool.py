#!/usr/bin/env python3
"""Pick the pool of non-osmotic Pfam families for the random-panel baseline.

The prespecified plan (README.md, "Plan provenance") compares the osmo panel against "a size-matched random
panel of 43 non-osmotic families". One random panel is one draw, so instead
this builds a POOL of families and analyze_extended.py draws many 43-family
panels from it, giving a null distribution for where the real panel falls.

Two subcommands, run in order:

  candidates  stream Pfam-A.hmm.gz, drop excluded families, draw a seeded
              random set of candidates, write candidates.hmm. Then search it
              against prokaryotic proteomes (hmmsearch --cut_ga --tblout).
  pool        read that tblout, keep candidates hit in >= --min-genomes of
              the screening proteomes, draw a seeded pool stratified by
              prevalence, write pool.hmm and pool_families.tsv.

Stratified pool: random Pfam families that occur in prokaryotes at all are
mostly rare (in the 2026-11 screen, ~75% hit under a quarter of genomes),
while the osmo panel is mostly mid-prevalence. A uniform draw would leave the
prevalence-matched null picking the same few mid-prevalence families every
time. So prevalence (share of screening proteomes hit) is binned in tenths,
and each bin gets --per-family x as many pool families as the osmo panel has
there (--osmo-matrix: a feature_matrix.tsv from the same screening genomes),
or all eligible ones if fewer exist.

Why screen: most of Pfam is eukaryotic or viral. A random draw from all of it
would fill panels with families absent from every genome here -- a "random
panel" that is mostly zeros says nothing about whether ANY gene content
predicts salt tolerance. Screening against the 43 pilot proteomes keeps
families that occur in BacDive-type prokaryotes at all; prevalence matching
in analyze_extended.py handles the rest.

Excluded: every Pfam accession in families.yaml (positives and negative
queries alike), plus any family whose NAME or DESC matches OSMO_TERMS.

Usage:
  python build_random_pool.py candidates --pfam Pfam-A.hmm.gz \\
      --families families.yaml --out-dir <dir> [--n 2000 --seed 1]
  hmmsearch --cut_ga --cpu 8 --tblout <dir>/screen.tbl -o /dev/null \\
      <dir>/candidates.hmm <pilot proteomes concatenated>
  python build_random_pool.py pool --pfam Pfam-A.hmm.gz --out-dir <dir> \\
      --screen <dir>/screen.tbl --osmo-matrix <pilot>/feature_matrix.tsv \\
      --n-genomes 43 [--min-genomes 2 --per-family 4 --seed 2]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import random
import re
from pathlib import Path

# Case-insensitive; matched against "NAME DESC". Deliberately broad: a family
# wrongly excluded costs nothing, one wrongly kept contaminates the null.
OSMO_TERMS = re.compile(
    r"osmo|betaine|choline|carnitine|ecto(ine)?\b|ectoine|trehalose|proline|"
    r"glycine.betaine|compatible.solute|mannosylglycerate|glucosylglycerol|"
    r"glucosylglycerate|sucrose|glutathione|potassium|sodium|\bK\+|\bNa\+|"
    r"cation|antiport|symport|mechanosensitive|\bMscS|\bMscL|\bTrk|\bKtr|"
    r"\bKdp|\bMrp|\bNha|\bBCCT|\bOpu|\bPro[PVWX]\b|\bBet[ABLT]\b|\bOts[AB]\b|"
    r"\bEct[ABCD]\b|cold.shock|\bCsp",
    re.IGNORECASE,
)


def read_hmms(path: Path):
    """Yield (acc_without_version, name, desc, record_text) per HMM."""
    with gzip.open(path, "rt") as fh:
        buf, acc, name, desc = [], None, None, ""
        for line in fh:
            buf.append(line)
            if line.startswith("NAME "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC "):
                acc = line.split(None, 1)[1].strip().split(".")[0]
            elif line.startswith("DESC "):
                desc = line.split(None, 1)[1].strip()
            elif line.startswith("//"):
                yield acc, name, desc, "".join(buf)
                buf, acc, name, desc = [], None, None, ""


def cmd_candidates(args) -> None:
    panel = set(re.findall(r"PF\d{5}", Path(args.families).read_text()))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    keep, n_total, n_panel, n_terms = [], 0, 0, 0
    for acc, name, desc, rec in read_hmms(Path(args.pfam)):
        n_total += 1
        if acc in panel:
            n_panel += 1
            continue
        if OSMO_TERMS.search(f"{name} {desc}"):
            n_terms += 1
            continue
        if "\nGA " not in rec:
            continue
        keep.append((acc, name, desc, rec))

    rng = random.Random(args.seed)
    chosen = rng.sample(keep, min(args.n, len(keep)))
    with open(out / "candidates.hmm", "w") as fh:
        for *_, rec in chosen:
            fh.write(rec)
    with open(out / "candidates.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pfam", "name", "desc"])
        w.writerows((a, n, d) for a, n, d, _ in chosen)
    print(f"Pfam-A families           : {n_total}")
    print(f"excluded, in families.yaml: {n_panel} (of {len(panel)} accessions listed)")
    print(f"excluded, osmo terms      : {n_terms}")
    print(f"eligible                  : {len(keep)}")
    print(f"candidates drawn (seed {args.seed}) : {len(chosen)} -> {out/'candidates.hmm'}")


def cmd_pool(args) -> None:
    out = Path(args.out_dir)
    cands = {r["name"]: r for r in csv.DictReader(open(out / "candidates.tsv"), delimiter="\t")}

    # tblout: target(protein) - query(model) ... ; protein ids from the pilot
    # proteomes are "<bacdive_id>|<prodigal id>" (see run_random_panel.sh)
    genomes_hit: dict[str, set[str]] = {}
    for line in open(args.screen):
        if line.startswith("#"):
            continue
        f = line.split()
        model, protein = f[2], f[0]
        genomes_hit.setdefault(model, set()).add(protein.split("|", 1)[0])

    eligible = sorted(m for m in cands if len(genomes_hit.get(m, ())) >= args.min_genomes)
    screened = {p.split("|", 1)[0] for line in open(args.screen) if not line.startswith("#")
                for p in [line.split()[0]]}
    n_genomes = args.n_genomes or len(screened)
    rng = random.Random(args.seed)

    def bin_of(prev: float) -> int:
        return min(int(prev * 10), 9)

    meta = {"bacdive_id", "species", "genus", "family", "domain", "strategy",
            "nacl_max_pct", "censored", "nogrowth_at_pct", "asm_level",
            "total_proteins", "n_families_present", "total_osmo_copies"}
    rows = list(csv.DictReader(open(args.osmo_matrix), delimiter="\t"))
    osmo = [c for c in rows[0] if c not in meta]
    want = [0] * 10
    for f in osmo:
        want[bin_of(sum(1 for r in rows if float(r[f]) > 0) / len(rows))] += 1

    by_bin: dict[int, list[str]] = {b: [] for b in range(10)}
    for m in eligible:
        by_bin[bin_of(len(genomes_hit[m]) / n_genomes)].append(m)
    pool = []
    print(f"{'prevalence':<12}{'osmo':>6}{'eligible':>10}{'pool':>6}")
    for b in range(10):
        k = min(len(by_bin[b]), want[b] * args.per_family)
        pool += rng.sample(by_bin[b], k)
        print(f"{b*10:>3}-{b*10+10:<3}%    {want[b]:>6}{len(by_bin[b]):>10}{k:>6}")
    pool = sorted(pool)

    wanted = set(pool)
    with open(out / "pool.hmm", "w") as fh:
        for acc, name, desc, rec in read_hmms(Path(args.pfam)):
            if name in wanted:
                fh.write(rec)
    with open(out / "pool_families.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pfam", "name", "desc", "screen_genomes_hit"])
        for m in pool:
            w.writerow([cands[m]["pfam"], m, cands[m]["desc"], len(genomes_hit.get(m, ()))])
    print(f"candidates                       : {len(cands)}")
    print(f"hit in >= {args.min_genomes} screening proteomes  : {len(eligible)}")
    print(f"pool drawn (seed {args.seed}, {args.per_family} per osmo family per bin): "
          f"{len(pool)} -> {out/'pool.hmm'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("candidates")
    c.add_argument("--pfam", required=True)
    c.add_argument("--families", required=True)
    c.add_argument("--out-dir", required=True)
    c.add_argument("--n", type=int, default=2000)
    c.add_argument("--seed", type=int, default=1)
    p = sub.add_parser("pool")
    p.add_argument("--pfam", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--screen", required=True)
    p.add_argument("--min-genomes", type=int, default=2)
    p.add_argument("--osmo-matrix", required=True,
                   help="feature_matrix.tsv for the screening genomes (the pilot's)")
    p.add_argument("--per-family", type=int, default=4)
    p.add_argument("--n-genomes", type=int, default=0,
                   help="screening proteome count; default: those with any hit")
    p.add_argument("--seed", type=int, default=2)
    args = ap.parse_args()
    {"candidates": cmd_candidates, "pool": cmd_pool}[args.cmd](args)


if __name__ == "__main__":
    main()
