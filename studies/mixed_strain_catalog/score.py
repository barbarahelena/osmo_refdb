#!/usr/bin/env python3
"""
score.py -- the pooled sample scored per genome group for DIAMOND alone and for
the catalog path under several settle rules.

Per mate: truth (truth.tsv: the annotated panel family it overlaps, or
background), DIAMOND's best hit from the run without a catalog (every mate),
and the catalog placements from the SAM (primary + XA). A mate's call under a
rule is the catalog's decision where the catalog settles it, else its DIAMOND
best hit, as compare_real_recall.py scores the catalog path (no report levels).
DIAMOND decides each mate on its own, so a fallback mate's call does not depend
on which other mates the catalog settled.

Methods:
  diamond_only          no catalog: every mate's DIAMOND best hit.
  current               osmotool 0.11.1's settle rule (settle_records): a mate
                        whose own alignment passes (>= 0.97 identity, >= 45 bp)
                        is settled; all passing placements are labelled
                        (annotated gene: its family, on a fused gene the family
                        whose envelope holds the read; anything else "drop");
                        one label -> call or drop, several -> ambiguous (no
                        call). A mate with a fall-through family among its
                        labels (catalog_reported "no") goes to DIAMOND.
  ambiguous_to_diamond  as current, but an ambiguous mate goes to DIAMOND.
  plus_fragments        as ambiguous_to_diamond, and a homolog gene whose bits
                        per covered HMM position reach the family's cutoff /
                        model length (a fragment of a true gene) counts as that
                        family, but only when every passing placement then
                        agrees on one family; it never makes a mate ambiguous.

Scored families: the release's panel minus the families profile does not
report (profile_excluded_families.txt); a call to any other label is no call.

Outputs in work/<tag>/: scores_by_group.tsv, scores_by_family.tsv,
decisions_by_group.tsv, and the current rule's decisions per mate
(settled_current.tsv) for checking against osmotool's own catalog_hits.tsv.
"""
from __future__ import annotations

import argparse, collections, csv, gzip, importlib.util, re
from pathlib import Path

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
ap = argparse.ArgumentParser()
for a in ("work", "release", "catalog", "tag", "compare"):
    ap.add_argument(f"--{a}", required=True)
args = ap.parse_args()
W, REL, CAT = Path(args.work), Path(args.release), Path(args.catalog)
OUT = W / args.tag

spec = importlib.util.spec_from_file_location("crr", args.compare)
crr = importlib.util.module_from_spec(spec); spec.loader.exec_module(crr)

groups = {r["genome"]: r["group"] for r in csv.DictReader(open(W.parent / "genomes.tsv"), delimiter="\t")}
truth, genome_of = {}, {}
for r in csv.DictReader(open(W / "truth.tsv"), delimiter="\t"):
    truth[r["read_id"]] = r["family"]; genome_of[r["read_id"]] = r["genome"]

excluded = {l.strip() for l in open(REL / "osmo_refdb.profile_excluded_families.txt") if l.strip()}
cut = {r["family"]: float(r["cutoff_bits"]) for r in csv.DictReader(open(REL / "hmms/cutoff_manifest.tsv"), delimiter="\t")}
panel = {f for f in cut if not f.endswith("_decoy")} - excluded
fallthrough = set()
with open(REL / "osmo_refdb.profile_report_levels.tsv") as fh:
    for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
        if r["catalog_reported"] == "no":
            fallthrough.update(r["members"].split(","))

# catalog annotation: label, status, envelopes (nt, 1-based on the gene)
ann = {}
for r in csv.DictReader(open(CAT / "catalog.annotation.tsv"), delimiter="\t"):
    env = (int(r["env_nt_start"]), int(r["env_nt_end"])) if r["env_nt_start"] else None
    others = []
    for item in filter(None, (r.get("other_envelopes") or "").split(",")):
        fam, span = item.rsplit(":", 1); a, b = span.split("-"); others.append((fam, int(a), int(b)))
    ann[r["gene"]] = (r["label"], r["status"], env, others, r["best_model"])

# fragments: homologs whose score density against their best model reaches cutoff / model length
cov, score, qlen = collections.defaultdict(int), {}, {}
with gzip.open(CAT / "catalog.hmmsearch.domtblout.gz", "rt") as fh:
    for l in fh:
        if l[0] == "#": continue
        p = l.split(); g, model = p[0], p[3]
        a = ann.get(g)
        if not a or a[1] != "homolog" or a[4] != model: continue
        cov[g] += int(p[16]) - int(p[15]) + 1; score[g] = float(p[7]); qlen[model] = int(p[5])
fragment = {g: ann[g][4] for g in cov if ann[g][4] in cut and score[g] / cov[g] >= cut[ann[g][4]] / qlen[ann[g][4]]}

def placement(gene, pos, cigar, nm):
    """(passes, label) as classify_placement; label "drop" for a non-call."""
    cols = al = span = 0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=X": cols += n; al += n; span += n
        elif op == "I": cols += n; al += n
        elif op in "DN": cols += n; span += n
    if not cols or 1 - nm / cols < 0.97 or al < 45:
        return False, None, None
    a = ann.get(gene)
    if a is None: return True, "drop", None
    label, status, env, others, _ = a
    if status != "annotated" or env is None:
        return True, "drop", fragment.get(gene)
    end = pos + span - 1
    ov = lambda x, y: max(0, min(end, y) - max(pos, x) + 1)
    own = ov(*env)
    for fam, x, y in others:
        o = ov(x, y)
        if o >= 0.5 * span and o > own:
            return True, fam, None
    return True, label, None

# per mate: catalog decisions under each rule; None = goes to DIAMOND
RULES = ("current", "ambiguous_to_diamond", "plus_fragments")
settled = {r: {} for r in RULES}       # rid -> ("call", fam) | ("drop",) | ("ambiguous",)
with open(OUT / "sample.sam") as fh:
    for line in fh:
        if line[0] == "@": continue
        f = line.rstrip("\n").split("\t"); flag = int(f[1])
        if flag & 0x900 or flag & 4: continue
        rid = f"{f[0]}/{1 if flag & 0x40 else 2}"
        nm = next(int(t[5:]) for t in f[11:] if t.startswith("NM:i:"))
        ok, lab, frag = placement(f[2], int(f[3]), f[5], nm)
        if not ok: continue
        pls = [(lab, frag)]
        for e in next((t[5:] for t in f[11:] if t.startswith("XA:Z:")), "").split(";"):
            p = e.split(",")
            if len(p) == 4:
                ok, lab, frag = placement(p[0], abs(int(p[1])), p[2], int(p[3]))
                if ok: pls.append((lab, frag))
        labels = {l for l, _ in pls}
        if labels & fallthrough:
            continue                                  # every rule: to DIAMOND
        if len(labels) > 1:
            settled["current"][rid] = ("ambiguous",)
        else:
            settled["current"][rid] = ("call", next(iter(labels))) if labels != {"drop"} else ("drop",)
        settled["ambiguous_to_diamond"][rid] = settled["current"][rid]
        loose = {fr if (l == "drop" and fr) else l for l, fr in pls}
        settled["plus_fragments"][rid] = (("call", next(iter(loose))) if len(loose) == 1 and loose != {"drop"}
                                  else settled["current"][rid])
        if loose & fallthrough:
            settled["plus_fragments"][rid] = settled["current"][rid]

diamond = crr.load_diamond_best_hits(W / "diamond" / "sample.blastx.tsv")
def calls(rule):
    if rule == "diamond_only":
        return {k: v[0] for k, v in diamond.items()}
    out = {}
    for rid in truth:
        d = settled[rule].get(rid)
        if d is None or (d[0] == "ambiguous" and rule != "current"):
            if rid in diamond: out[rid] = diamond[rid][0]
        elif d[0] == "call":
            out[rid] = d[1]
    return out

methods = ("diamond_only",) + RULES
called = {m: {k: v for k, v in calls(m).items() if v in panel} for m in methods}

with open(OUT / "settled_current.tsv", "w") as fh:
    fh.write("read_id\tdecision\tlabel\n")
    for rid, d in settled["current"].items():
        fh.write(f"{rid}\t{d[0]}\t{d[1] if len(d) > 1 else ''}\n")

order = ["in", "partial", "out_related", "out_distant", "all"]
def grp(rid): return groups[genome_of[rid]]
rows = []
for m in methods:
    c = called[m]
    for g in order:
        ids = [r for r in truth if g == "all" or grp(r) == g]
        pos = [r for r in ids if truth[r] in panel]
        tp = sum(c.get(r) == truth[r] for r in pos)
        wrong_fam = sum(r in c and c[r] != truth[r] for r in pos)
        fp_bg = sum(r in c for r in ids if truth[r] not in panel)
        fp = wrong_fam + fp_bg
        rows.append({"method": m, "group": g, "panel_mates": len(pos), "tp": tp, "fn": len(pos) - tp,
                     "fp": fp, "fp_wrong_family": wrong_fam, "fp_background": fp_bg,
                     "recall": round(tp / len(pos), 3) if pos else "",
                     "precision": round(tp / (tp + fp), 3) if tp + fp else "",
                     "settled_by_catalog": "" if m == "diamond_only" else sum(r in settled[m] for r in ids)})
with open(OUT / "scores_by_group.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

# sample-level: reads per family, called vs true
fams = sorted(panel)
with open(OUT / "scores_by_family.tsv", "w") as fh:
    fh.write("family\ttrue_mates\t" + "\t".join(f"{m}_called\t{m}_tp" for m in methods) + "\n")
    for fam in fams:
        t = sum(v == fam for v in truth.values())
        cells = []
        for m in methods:
            c = called[m]; n = sum(v == fam for v in c.values()); tp = sum(v == fam and truth[r] == fam for r, v in c.items())
            cells += [str(n), str(tp)]
        if t or any(int(x) for x in cells):
            fh.write(f"{fam}\t{t}\t" + "\t".join(cells) + "\n")

# catalog decisions per group (current rule) and for panel mates only
dec = collections.Counter()
for rid in truth:
    d = settled["current"].get(rid, ("to_DIAMOND",))[0]
    dec[(grp(rid), "panel" if truth[rid] in panel else "background", d)] += 1
with open(OUT / "decisions_by_group.tsv", "w") as fh:
    fh.write("group\tmates\tdecision\tn\n")
    for (g, kind, d), n in sorted(dec.items()):
        fh.write(f"{g}\t{kind}\t{d}\t{n}\n")

# presentation table: groups as rows, methods as columns, recall and false calls side by side
NAMES = {"diamond_only": "DIAMOND only (no catalogue)", "current": "Catalogue, current rule",
         "ambiguous_to_diamond": "Catalogue, ambiguous reads to DIAMOND",
         "plus_fragments": "Catalogue, ambiguous to DIAMOND + fragment rule"}
n_genomes = collections.Counter(groups.values())
GROUPS = {"in": f"Strains in the catalogue ({n_genomes['in']} genomes)",
          "partial": f"Strains partly in the catalogue ({n_genomes['partial']})",
          "out_related": f"Not in the catalogue, same genus as a catalogued strain ({n_genomes['out_related']})",
          "out_distant": f"Not in the catalogue, distant ({n_genomes['out_distant']})",
          "all": f"Whole sample ({len(groups)})"}
cell = {(r["method"], r["group"]): r for r in rows}
with open(OUT / "results_table.md", "w") as fh:
    for title, key, fmt in (("Reads from panel genes that get the right family (recall)", "recall",
                             lambda r: f"{r['recall']:.2f} ({r['tp']}/{r['panel_mates']})"),
                            ("False calls (reads given a family they are not from)", "fp",
                             lambda r: f"{r['fp']:,}")):
        fh.write(f"**{title}**\n\n| Strains | " + " | ".join(NAMES[m] for m in methods) + " |\n")
        fh.write("|---|" + "---:|" * len(methods) + "\n")
        for g in order:
            fh.write(f"| {GROUPS[g]} | " + " | ".join(fmt(cell[(m, g)]) for m in methods) + " |\n")
        fh.write("\n")

print(f"{len(fragment)} fragment homologs; fall-through families: {','.join(sorted(fallthrough))}")
print(f"{'method':14s} {'group':12s} {'panel':>6s} {'tp':>5s} {'fp':>6s} {'(wrong fam / bg)':>17s} {'recall':>7s} {'prec':>6s} {'settled':>8s}")
for r in rows:
    print(f"{r['method']:14s} {r['group']:12s} {r['panel_mates']:6d} {r['tp']:5d} {r['fp']:6d} "
          f"{r['fp_wrong_family']:>8d} / {r['fp_background']:<6d} {r['recall']!s:>7s} {r['precision']!s:>6s} {r['settled_by_catalog']!s:>8s}")
