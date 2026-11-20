#!/usr/bin/env python3
"""
05b_build_hmm_decoys.py — build HMM decoy models: profiles of a confuser
class that compete in osmotool `annotate`'s best-hit contest and are then
dropped from its output, the HMM counterpart of 08a_build_decoy_refs.py's
DIAMOND decoys.

Motivating case: opuBA. The v9 BacDive run called opuBA in 64% of genomes
and 48% of salt-in strains. Most of those proteins were peptide ABC
ATPases (OppD/OppF/DppD-type, C-terminal oligo_HPY domain PF08352), which
no opuAA/BA/CA positive carries. Adding them to the opu ATPase negative
pool did not help opuBA: they overlap its 23 calibration positives in
score, so the F1 sweep kept the GA at ~295 bits. A model of the peptide
ATPases themselves scores those proteins far higher than any opu model
does, so they go to the decoy and are not counted; in a 484-genome
prototype it took all 744 such proteins and none of the 128 held-out opu
positives (each scored >= 57 bits higher on its own model).

How it plugs in, with no new scoring logic anywhere:
  * annotate: hmmsearch --cut_ga over the pressed db, best hit per protein,
    then labels in osmo_refdb.annotate_excluded_families.txt are dropped
    (08c lists every hmm_decoys name there and in the profile list).
  * benchmark: 10_run_benchmark.sh scans reads against the same pressed db;
    11_compute_metrics.py counts a best hit to a *_decoy label as NO_CALL.
  * 07_press_hmms.sh concatenates hmms/*.hmm, so the decoy is included.
  * 06_calibrate_cutoffs.py only iterates families.yaml `families`, so a
    decoy's GA is the fixed `ga` given here. It only has to be below the
    competing families' cutoffs: a decoy can only win a protein it scores
    higher than they do.

Leakage: any accession already present in another refs/*.faa file (every
family's positives and negatives, train or test, and the random
background) is dropped before building, so no benchmark read comes from a
sequence the decoy was trained on.

`from_negatives_of: <family>` builds the decoy from that family's
negative.train instead of a fetch: for a family whose named confusers are
already its whole hard-negative pool (betL: every UniProt betT/caiT is in
it), so no independent set exists. That is the HMM counterpart of
decoy_from_negatives' DIAMOND decoy, on the same split: negative.test,
which the benchmark scores, is disjoint from it. Any accession also present
in a test split is dropped anyway.

`decoy_aware_ga: true` makes the competing families' GA calibration
decoy-aware (06_calibrate_cutoffs.py): their positive.calib and
negative.train are scored against the decoy here
(hmms/scores/<family>.<decoy>.{positive,negative}.tblout), and a sequence
the decoy outscores is not called as the family, as in annotate's best-hit
contest. Caveat for a from_negatives_of decoy: it was trained on that very
negative.train, so it wins those negatives more easily than unseen ones;
the cutoff it allows is therefore optimistic, and only the test split
(negative.test) shows what it costs.

Usage:
  python 05b_build_hmm_decoys.py --refs refs --aln alignments --hmms hmms \\
      --families families.yaml [--threads 4]
Input:  refs/<decoy>.hmm_decoy.raw.faa   (01_fetch_refs.py)
Output: refs/<decoy>.hmm_decoy.faa, refs/<decoy>.hmm_decoy.clustered.faa,
        <aln>/<decoy>.trimmed.aln, <hmms>/<decoy>.hmm (with GA line)
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent


def load_set_ga_line():
    """06_calibrate_cutoffs.py isn't importable by name (leading digit)."""
    spec = importlib.util.spec_from_file_location("calibrate", HERE / "06_calibrate_cutoffs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.set_ga_line


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


def accession(header: str) -> str:
    """'<tag>|UniProtID|Organism' -> 'UniProtID'."""
    parts = header.split()[0].split("|")
    return parts[1] if len(parts) > 1 else parts[0]


def test_accessions(refs: Path) -> set[str]:
    """Accessions the benchmark scores: every family's positive.test and
    negative.test, and the random background."""
    out: set[str] = set()
    for pattern in ("*.positive.test.faa", "*.negative.test.faa", "random_background.faa"):
        for path in refs.glob(pattern):
            out.update(accession(h) for h, _ in read_fasta(path))
    return out


def used_accessions(refs: Path) -> set[str]:
    used: set[str] = set()
    for path in refs.glob("*.faa"):
        if ".hmm_decoy." in path.name:
            continue
        used.update(accession(h) for h, _ in read_fasta(path))
    return used


def run(cmd: list[str], stdout: Path | None = None) -> None:
    if stdout is None:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    else:
        with open(stdout, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.DEVNULL)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", type=Path, default=Path("refs"))
    ap.add_argument("--aln", type=Path, default=Path("alignments"))
    ap.add_argument("--hmms", type=Path, default=Path("hmms"))
    ap.add_argument("--families", type=Path, default=Path("families.yaml"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--reuse-alignments", action="store_true",
                    help="Keep an existing decoy clustering + alignment when the "
                         "decoy's input sequences are unchanged.")
    args = ap.parse_args()

    with open(args.families) as fh:
        decoys = yaml.safe_load(fh).get("hmm_decoys") or []
    if not decoys:
        print("No hmm_decoys in families.yaml; nothing to build.")
        return

    set_ga_line = load_set_ga_line()
    used = used_accessions(args.refs)
    tested = test_accessions(args.refs)
    args.aln.mkdir(parents=True, exist_ok=True)
    args.hmms.mkdir(parents=True, exist_ok=True)

    for decoy in decoys:
        name = decoy["name"]
        if not name.endswith("_decoy"):
            # 11_compute_metrics.py recognises decoy calls by this suffix.
            raise SystemExit(f"ERROR: hmm_decoys name {name!r} must end in '_decoy'")
        source_family = decoy.get("from_negatives_of")
        if source_family:
            src = args.refs / f"{source_family}.negative.train.faa"
            records = read_fasta(src)
            kept = [(h, s) for h, s in records if accession(h) not in tested]
            print(f"[{name}] {len(records)} from {src.name}, {len(records) - len(kept)} dropped "
                  f"(accession in a test split), {len(kept)} kept")
        else:
            raw = args.refs / f"{name}.hmm_decoy.raw.faa"
            records = read_fasta(raw)
            kept = [(h, s) for h, s in records if accession(h) not in used]
            print(f"[{name}] {len(records)} fetched, {len(records) - len(kept)} dropped "
                  f"(accession already in another refs set), {len(kept)} kept")
        if len(kept) < 50:
            raise SystemExit(f"ERROR: [{name}] only {len(kept)} sequences left to build from")

        filtered = args.refs / f"{name}.hmm_decoy.faa"
        clustered = args.refs / f"{name}.hmm_decoy.clustered.faa"
        raw_aln = args.aln / f"{name}.aln"
        trimmed = args.aln / f"{name}.trimmed.aln"
        new_filtered = "".join(f">{h}\n{s}\n" for h, s in kept)
        if (args.reuse_alignments and trimmed.exists() and clustered.exists()
                and filtered.exists() and filtered.read_text() == new_filtered):
            # Same input sequences as the build these files came from:
            # keep its clustering and alignment (multi-threaded MAFFT is
            # not reproducible run to run).
            print(f"[{name}] input unchanged, reusing {trimmed}")
        else:
            filtered.write_text(new_filtered)
            run(["cd-hit", "-i", str(filtered), "-o", str(clustered), "-c", "0.90", "-n", "5",
                 "-T", str(args.threads), "-M", "16000", "-d", "0"])
            run(["mafft", "--auto", "--thread", str(args.threads), str(clustered)], stdout=raw_aln)
            run(["trimal", "-in", str(raw_aln), "-out", str(trimmed), "-gt", str(decoy.get("trim_gt", 0.8))])
        n_clu = sum(1 for h, _ in read_fasta(clustered))

        hmm = args.hmms / f"{name}.hmm"
        run(["hmmbuild", "--amino", "-n", name, "--cpu", str(args.threads), str(hmm), str(trimmed)])
        set_ga_line(hmm, float(decoy["ga"]))
        leng = next(l.split()[1] for l in hmm.read_text().splitlines() if l.startswith("LENG"))
        print(f"[{name}] {n_clu} clusters -> {hmm} (LENG {leng}, GA {float(decoy['ga']):.1f}; "
              f"competes with {', '.join(decoy.get('competes_with', [])) or 'unspecified'})")

        if decoy.get("decoy_aware_ga"):
            # 06_calibrate_cutoffs.py: the competing families' calibration
            # sequences, scored against this decoy
            score_dir = args.hmms / "scores"
            score_dir.mkdir(exist_ok=True)
            for fam in decoy.get("competes_with", []):
                for kind, split in (("positive", "positive.calib"), ("negative", "negative.train")):
                    faa = args.refs / f"{fam}.{split}.faa"
                    if not faa.exists() or faa.stat().st_size == 0:
                        continue
                    out = score_dir / f"{fam}.{name}.{kind}.tblout"
                    run(["hmmsearch", "--cpu", str(args.threads), "--tblout", str(out),
                         "--noali", str(hmm), str(faa)])
                print(f"[{name}] scored {fam}'s calibration sets for decoy-aware GA")


if __name__ == "__main__":
    main()
