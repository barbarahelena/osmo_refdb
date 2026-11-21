#!/usr/bin/env python3
"""
package_catalog.py -- build a gene catalog into a release as
catalogs/<name>/, for `osmotool profile --catalog <name>`.

The catalog's annotation is one unit with the release it is built for: it
comes from this release's HMMs and gathering cutoffs, so rebuild it for
every release (it changes whenever an HMM or cutoff does). It ships with
the index osmotool build-catalog --index neighbours makes: the genes with an
HMM hit (sub-cutoff homologs and decoys included, so their reads are
recognised and dropped) plus the genes without a hit that are close to one,
so reads from those competitors map to their own gene and are dropped too.
Small enough to ship, and on real gut genes as free of catalog false calls as
an index of the whole catalog (docs/validation_reports/catalog_index.md).

The work is done by `osmotool build-catalog` (translate, hmmsearch with the
release's HMMs, assign each gene as `annotate --method hmm` does, index the
whole catalog with bwa-mem3). This script adds what makes it a release
component: provenance (source, input checksums, release), a compressed copy
of the raw hmmsearch evidence, and SHA256SUMS.

Usage:
  python bin/package_catalog.py --release releases/v10 --name gmgc_gut \\
      --genes GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz \\
      --source "GMGC v1.0, human gut, prevalence > 0.5%" \\
      --source-url https://gmgc.embl.de/ --threads 32

Writes releases/<release>/catalogs/<name>/:
  catalog.fna.*                    bwa-mem3 index of the genes with an HMM hit
                                   and their neighbours (--index neighbours)
  catalog.annotation.tsv           every gene with an HMM hit: family (or
                                   homolog / decoy), HMM envelope in nt
  catalog.json                     build parameters, counts per family,
                                   provenance
  catalog.hmmsearch.domtblout.gz   raw hmmsearch domain table behind the
                                   annotation
  SHA256SUMS                       checksums of the files above
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", type=Path, default=None,
                    help="release directory (needed unless --built, --out and --release-name are given)")
    ap.add_argument("--release-name", default=None, help="default: the --release directory's name")
    ap.add_argument("--name", required=True,
                    help="catalog name, used as `profile --catalog NAME` (e.g. gmgc_gut)")
    ap.add_argument("--genes", type=Path, required=True,
                    help="nucleotide CDS FASTA (gz ok), one in-frame gene per record")
    ap.add_argument("--proteins", type=Path, default=None,
                    help="protein FASTA with the same IDs (otherwise genes are translated)")
    ap.add_argument("--source", default="", help="what the catalog is, in words")
    ap.add_argument("--source-url", default="", help="where it was downloaded from")
    ap.add_argument("--index", choices=["hits", "all", "none"], default="hits",
                    help="hits (default): index only genes with an HMM hit, small enough "
                         "to ship with the release; all: every gene (large); none: "
                         "annotation only")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--osmotool", default="osmotool",
                    help="osmotool command (needs build-catalog)")
    ap.add_argument("--built", type=Path, default=None,
                    help="a catalog directory osmotool build-catalog has already written: "
                         "package that instead of running osmotool (the Nextflow pipeline does)")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the packaged catalog here instead of <release>/catalogs/<name>/ "
                         "(the Nextflow pipeline publishes it into the release)")
    ap.add_argument("--force", action="store_true",
                    help="replace an existing catalogs/<name>/")
    args = ap.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.name):
        sys.exit(f"--name {args.name!r}: use letters, digits, '_', '-', '.'")
    if args.release is None and not (args.built and args.out and args.release_name):
        ap.error("--release is required unless --built, --out and --release-name are given")
    release_files = [] if args.built else [args.release / "osmo_refdb.dmnd",
                                           args.release / "hmms" / "osmo_refdb.hmm"]
    for f in release_files + [args.genes] + ([args.proteins] if args.proteins else []):
        if not f.exists():
            sys.exit(f"ERROR: {f} not found")
    out = args.out or args.release / "catalogs" / args.name
    if out.exists():
        if not args.force:
            sys.exit(f"ERROR: {out} exists; pass --force to rebuild it")
        shutil.rmtree(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.built:
        if not (args.built / "catalog.json").exists():
            sys.exit(f"ERROR: {args.built}/catalog.json not found: not a built catalog")
        shutil.copytree(args.built, out)   # follows the link a staged input is
    else:
        cmd = [*args.osmotool.split(), "build-catalog", str(args.release), "--genes", str(args.genes),
               "-o", str(out), "--index", args.index, "--threads", str(args.threads)]
        if args.proteins:
            cmd += ["--proteins", str(args.proteins)]
        print("Running:", " ".join(cmd), flush=True)
        if subprocess.run(cmd).returncode != 0:
            sys.exit("ERROR: osmotool build-catalog failed")

    # raw evidence, compressed
    domtbl = out / "catalog.hmmsearch.domtblout"
    if domtbl.exists():
        with open(domtbl, "rb") as fin, gzip.open(str(domtbl) + ".gz", "wb") as fout:
            shutil.copyfileobj(fin, fout)
        domtbl.unlink()

    # provenance: the build machine's paths mean nothing in a release
    meta_path = out / "catalog.json"
    meta = json.loads(meta_path.read_text())
    release_name = args.release_name or args.release.name
    meta.update({
        "name": args.name,
        # build-catalog wrote the name of the directory it was given, which under Nextflow is the
        # staged copy ("db"); the release's own name replaces it
        "release": release_name,
        "osmo_refdb_release": release_name,
        "packaged": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": args.source,
        "source_url": args.source_url,
        "genes": args.genes.name,
        "genes_sha256": meta.get("genes_sha256") or sha256_of(args.genes),
        "proteins": args.proteins.name if args.proteins else "",
        "proteins_sha256": sha256_of(args.proteins) if args.proteins else "",
        "hmm_db": "hmms/osmo_refdb.hmm",
    })
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    files = sorted(p for p in out.iterdir() if p.is_file() and p.name != "SHA256SUMS")
    (out / "SHA256SUMS").write_text("".join(f"{sha256_of(p)}  {p.name}\n" for p in files))

    size = sum(p.stat().st_size for p in files)
    per_label = meta.get("genes_per_label", {})
    n_annot = sum(v for k, v in per_label.items() if k not in ("homolog", "decoy"))
    print(f"Catalog {args.name} packaged into {out}: {meta.get('n_genes_indexed')} genes indexed, "
          f"{n_annot} annotated to a family, {size / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
