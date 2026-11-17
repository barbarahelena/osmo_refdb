#!/usr/bin/env python3
"""
strip_ref_tag.py -- remove every sequence with a given source tag from a
refs/ directory copied from another release (REUSE_REFS_FROM), offline, and
re-cluster the families it touched.

Why: step 1d merged local study sequences (tag "<family>_study", from
extra_sequences/, which is gitignored and cannot be regenerated) into four
families' positives. It was retired (2026-11-17), but releases built before
still carry those sequences in their refs/. A rebuild from such refs strips
them here before the train/calib/test split.

Step 02 overwrites <family>.positive.faa with the clustered set, so a tagged
sequence that was a cluster representative may have absorbed other
sequences. Deleting it would silently drop them too. Instead this script:
  1. drops the tagged records from <family>.positive.faa;
  2. recovers the untagged members its clusters had absorbed, from
     <family>.positive.pre_length_filter.faa (UniProt entries before the
     length filter; step 02's input minus the extras). Members from other
     sources (RefSeq, "_refseq" tag) are not saved anywhere and are counted
     as lost, not recovered;
  3. clusters what is left again with CD-HIT at 90% (as step 02), and
     rewrites <family>.positive.faa, .positive.clustered.faa and its .clstr.
Writes refs/removed_tag_manifest.tsv (what was removed, recovered and lost,
per family), so nothing disappears unannounced.

Usage:
  python bin/strip_ref_tag.py --refs releases/v10-rebuild/refs --tag _study
Needs cd-hit on PATH. A refs/ directory without the tag is left untouched.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, header, seq = [], None, []
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


def token(header: str) -> str:
    return header.split()[0]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def clusters_of(clstr: Path) -> list[list[tuple[str, bool]]]:
    """[(member name, is representative)] per cluster, from a CD-HIT .clstr."""
    clusters: list[list[tuple[str, bool]]] = []
    for line in clstr.read_text().splitlines():
        if line.startswith(">Cluster"):
            clusters.append([])
        elif line.strip():
            name = re.search(r">(\S+?)\.\.\.", line).group(1)
            clusters[-1].append((name, line.rstrip().endswith("*")))
    return clusters


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, required=True)
    ap.add_argument("--tag", default="_study", help="source tag to remove, as in '<family>_study|'")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    marker = f"{args.tag}|"

    manifest = []
    for pos in sorted(args.refs.glob("*.positive.faa")):
        family = pos.name[:-len(".positive.faa")]
        records = read_fasta(pos)
        tagged = [h for h, _ in records if marker in h]
        if not tagged:
            continue
        kept = [(h, s) for h, s in records if marker not in h]

        absorbed: list[str] = []
        clstr = args.refs / f"{family}.positive.clustered.faa.clstr"
        if clstr.exists():
            for cluster in clusters_of(clstr):
                rep = next((n for n, is_rep in cluster if is_rep), None)
                if rep is not None and marker in rep:
                    absorbed += [n for n, is_rep in cluster if not is_rep and marker not in n]
        pre = args.refs / f"{family}.positive.pre_length_filter.faa"
        source = {token(h): (h, s) for h, s in read_fasta(pre)} if pre.exists() else {}
        recovered = [source[n] for n in absorbed if n in source]
        lost = Counter(re.search(r"_([A-Za-z]+)\|", n).group(1) if re.search(r"_([A-Za-z]+)\|", n)
                       else "other" for n in absorbed if n not in source)

        pool = {token(h): (h, s) for h, s in kept + recovered}
        with tempfile.TemporaryDirectory() as tmp:
            infa, out = Path(tmp) / "pool.faa", Path(tmp) / "clustered.faa"
            write_fasta(infa, list(pool.values()))
            if len(pool) >= 5:
                subprocess.run(["cd-hit", "-i", str(infa), "-o", str(out), "-c", "0.90", "-n", "5",
                                "-T", str(args.threads), "-M", "16000", "-d", "0"],
                               check=True, stdout=subprocess.DEVNULL)
                final = read_fasta(out)
                clustered_text, clstr_text = out.read_text(), Path(str(out) + ".clstr").read_text()
            else:
                final = list(pool.values())
                clustered_text, clstr_text = infa.read_text(), ""
        write_fasta(pos, final)
        (args.refs / f"{family}.positive.clustered.faa").write_text(clustered_text)
        clstr.write_text(clstr_text)
        for stale in ("positive.train.faa", "positive.calib.faa", "positive.test.faa"):
            (args.refs / f"{family}.{stale}").unlink(missing_ok=True)

        row = {"family": family, "tag": args.tag, "n_positives_before": len(records),
               "n_removed": len(tagged), "n_recovered_untagged": len(recovered),
               "n_lost_absorbed": sum(lost.values()),
               "lost_by_source": ",".join(f"{k}:{v}" for k, v in sorted(lost.items())),
               "n_positives_after": len(final)}
        manifest.append(row)
        print(f"[{family}] removed {row['n_removed']} {args.tag} sequences, recovered "
              f"{row['n_recovered_untagged']} they had absorbed, lost {row['n_lost_absorbed']} "
              f"({row['lost_by_source'] or '-'}); {row['n_positives_before']} -> {row['n_positives_after']}")

    if manifest:
        path = args.refs / "removed_tag_manifest.tsv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(manifest[0]), delimiter="\t")
            w.writeheader()
            w.writerows(manifest)
        print(f"Written: {path}")
    else:
        print(f"No sequences tagged {args.tag} in {args.refs}: nothing to strip.")


if __name__ == "__main__":
    main()
