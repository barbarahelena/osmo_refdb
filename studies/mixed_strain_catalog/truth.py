#!/usr/bin/env python3
"""
truth.py -- per-mate truth for the pooled sample: the family of the annotated
panel gene the mate overlaps by >= 20 bp (real_read_truth.py's rule), else
"background", plus the genome the mate came from.

wgsim names are <contig>_<start>_<end>_<e1>_<e2>_<n>, with start/end the
fragment's outer coordinates (1-based). Mate 1 is at one end of the fragment
and mate 2 at the other, but which end depends on the strand wgsim drew, which
the name does not say; it is read from the FASTQ by matching the mate's
sequence against both ends of the fragment (the end with fewer mismatches).

Writes work/truth.tsv: read_id (<name>/<mate>), genome, family.
"""
import csv, sys
from pathlib import Path

MIN_OVERLAP = 20
W = Path(sys.argv[1])
COMP = str.maketrans("ACGTN", "TGCAN")

def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name: seqs[name] = "".join(buf)
    return seqs

genomes = [r["genome"] for r in csv.DictReader(open(W.parent / "genomes.tsv"), delimiter="\t")]
genes = {}   # contig -> [(start, end, family)]
for g in genomes:
    for r in csv.DictReader(open(W / "annotate" / f"{g}.gene_coordinates.tsv"), delimiter="\t"):
        genes.setdefault(r["contig"], []).append((int(r["start"]), int(r["end"]), r["gene"]))

def mismatches(a, b):
    return sum(x != y for x, y in zip(a, b))

n_total = n_panel = 0
with open(W / "truth.tsv", "w") as out:
    out.write("read_id\tgenome\tfamily\n")
    for g in genomes:
        seqs = fasta(W / "genomes" / f"{g}.fna")
        for mate in (1, 2):
            with open(W / "reads" / f"{g}.R{mate}.fq") as fq:
                while True:
                    head = fq.readline()
                    if not head: break
                    seq = fq.readline().strip(); fq.readline(); fq.readline()
                    name = head[1:].split()[0].rsplit("/", 1)[0]
                    contig, start, end = name.rsplit("_", 5)[:3]
                    start, end, L = int(start), int(end), len(seq)
                    ref = seqs[contig]
                    left = ref[start - 1:start - 1 + L]
                    right = ref[end - L:end].translate(COMP)[::-1]
                    lo = start if mismatches(seq, left) <= mismatches(seq, right) else end - L + 1
                    hi = lo + L - 1
                    fam = "background"
                    for a, b, f in genes.get(contig, ()):
                        if min(hi, b) - max(lo, a) + 1 >= MIN_OVERLAP:
                            fam = f; break
                    out.write(f"{name}/{mate}\t{g}\t{fam}\n")
                    n_total += 1; n_panel += fam != "background"
print(f"{n_total} mates, {n_panel} overlap a panel gene")
