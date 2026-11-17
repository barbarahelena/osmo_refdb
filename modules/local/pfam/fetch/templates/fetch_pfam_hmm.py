#!/usr/bin/env python3
"""
fetch_pfam_hmm.py -- a family's curated Pfam model (families.yaml `pfam_model`) from InterPro, used
in place of an hmmbuild model. Keeps Pfam's own GA line; NAME and ACC are rewritten to the family's
name: hmmpress's index needs both unique across the combined database, and several families can adopt
the same Pfam accession (trkH/ktrB/ktrD -> PF02386).

Usage: fetch_pfam_hmm.py <accession> <family> <out.hmm>
"""
import gzip
import re
import sys

import requests


def main() -> None:
    accession, family, hmm_file = sys.argv[1], sys.argv[2], sys.argv[3]
    r = requests.get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/",
                     params={"annotation": "hmm"}, timeout=30)
    r.raise_for_status()
    text = gzip.decompress(r.content).decode()
    text = re.sub(r"^NAME  .*$", f"NAME  {family}", text, count=1, flags=re.MULTILINE)
    if re.search(r"^ACC   ", text, flags=re.MULTILINE):
        text = re.sub(r"^ACC   .*$", f"ACC   {family}", text, count=1, flags=re.MULTILINE)
    with open(hmm_file, "w") as fh:
        fh.write(text)
    print(f"[{family}] curated Pfam model {accession} -> {hmm_file}")


if __name__ == "__main__":
    main()
