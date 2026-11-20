#!/usr/bin/env python3
"""
write_provenance.py -- record what produced a release, next to the release.

Writes <release>/PROVENANCE.json:
  * code: the osmo_refdb git commit and branch, and whether any files had
    been changed since that commit without being committed (and which), so a
    build is either tied to one exact version of the code or says it is not;
    osmotool's version, git commit (from pip's direct_url.json) and the path it
    was imported from, which shows a PYTHONPATH working tree shadowing the
    installed one;
  * environment: machine architecture and CPU model, host, container runtime,
    the image's SHA-256 (SIF_SHA256, set by the sbatch script), the hash of the
    explicit conda list baked into the image, the version of every tool the
    pipeline calls;
  * inputs: families.yaml checksum, the UniProt release(s) and fetch dates in
    refs/manifest.tsv, the release refs/ were reused from (REUSE_REFS_FROM),
    pipeline parameters from the environment (SPLIT_MODE, REUSE_ALIGNMENTS,
    ...), and a checksum over refs/ and alignments/;
  * outputs: checksums of every HMM (without the DATE and COM lines HMMER
    writes), the DIAMOND database, cutoffs, cascade config and report levels.

This is a record, not a guarantee: it says what ran so a rebuild can be
compared (hpc/verify_reproducible.py) and a difference traced.

Usage:
  python pipeline/write_provenance.py --release releases/v10 [--families families.yaml]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

TOOLS = {
    "mafft": (["mafft", "--version"], r"v?\d[\w.]*\s*\(.*?\)|v?\d[\w.]*"),
    "trimal": (["trimal", "--version"], r"trimAl.*"),
    "hmmbuild": (["hmmbuild", "-h"], r"HMMER\s+[\d.]+[^;]*"),
    "diamond": (["diamond", "--version"], r"diamond version [\d.]+"),
    "cd-hit": (["cd-hit", "-h"], r"CD-HIT version [\d.]+"),
    "prodigal": (["prodigal", "-v"], r"Prodigal V[\d.]+"),
    "orfm": (["orfm", "--version"], r"[\d.]+"),
    "wgsim": (["wgsim"], r"Version: [\d.]+\S*"),
    "iss": (["iss", "--version"], r"[\d.]+"),
    "bwa-mem3": (["bwa-mem3", "version"], r"\d+\.\d+[\w.\-]*"),
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def hmm_sha256(path: Path) -> str:
    """Checksum of an .hmm file without the DATE and COM lines HMMER writes
    (a timestamp and the hmmbuild command line, with file paths): two builds
    of the same alignment differ in exactly those two lines."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for line in fh:
            if not line.startswith((b"DATE", b"COM")):
                h.update(line)
    return h.hexdigest()


def tree_hash(directory: Path, pattern: str = "*") -> dict:
    """One checksum over the files of a directory: name + content hash of each."""
    files = sorted(p for p in directory.glob(pattern) if p.is_file())
    h = hashlib.sha256()
    for p in files:
        h.update(p.name.encode() + b"\0" + sha256_of(p).encode() + b"\n")
    return {"n_files": len(files), "sha256": h.hexdigest()} if files else {"n_files": 0}


def run(cmd: list[str], cwd: Path | None = None, strip: bool = True) -> str:
    """Output of a command. strip=False keeps leading whitespace, which
    `git status --porcelain` needs: its first column can be a space."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, timeout=60)
        text = r.stdout + r.stderr
        return text.strip() if strip else text.rstrip("\n")
    except (OSError, subprocess.SubprocessError):
        return ""


def git_info(repo: Path) -> dict:
    g = ["git", "-c", "safe.directory=*", "-C", str(repo)]
    commit = run(g + ["rev-parse", "HEAD"])
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        return {"commit": "unknown (not a git checkout)"}
    changed = [l[3:] for l in run(g + ["status", "--porcelain"], strip=False).splitlines()
               if not l[3:].startswith(("releases/", "cache/", ".claude/"))
               and not l[3:].endswith(".DS_Store")]
    return {"commit": commit, "branch": run(g + ["rev-parse", "--abbrev-ref", "HEAD"]),
            "describe": run(g + ["describe", "--tags", "--always"]),
            "uncommitted_changes": bool(changed), "uncommitted_files": changed[:50]}


def osmotool_info() -> dict:
    """The osmotool that is actually imported, not the one pip lists: a
    PYTHONPATH working tree shadows the installed package, and its metadata
    (version, direct_url) then describes the wrong thing. If the imported
    code lies outside site-packages, its own git checkout is recorded."""
    out: dict = {}
    try:
        import sysconfig
        import osmotool
        pkg = Path(osmotool.__file__).resolve().parent
        out["version"] = getattr(osmotool, "__version__", "unknown")
        out["imported_from"] = str(pkg)
        site = {Path(sysconfig.get_paths()[k]).resolve() for k in ("purelib", "platlib")}
        installed = any(s in pkg.parents for s in site)
        out["shadowed_by_pythonpath"] = not installed
        if installed:
            from importlib import metadata
            for f in metadata.distribution("osmotool").files or []:
                if f.name == "direct_url.json":
                    url = json.loads(f.read_text())
                    out["installed_from"] = url.get("url", "")
                    out["git_commit"] = url.get("vcs_info", {}).get("commit_id", "")
        else:
            out["working_tree_git"] = git_info(pkg.parent.parent)
    except Exception as exc:       # noqa: BLE001 -- a record must not crash a build
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def tool_versions() -> dict:
    versions = {}
    for name, (cmd, pattern) in TOOLS.items():
        text = run(cmd)
        m = re.search(pattern, text)
        versions[name] = m.group(0).strip() if m else ("not found" if not text else text.splitlines()[0][:80])
    return versions


def cpu_model() -> str:
    try:
        for line in open("/proc/cpuinfo"):
            if line.lower().startswith(("model name", "hardware", "cpu model")):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--families", type=Path, default=REPO / "families.yaml")
    ap.add_argument("--repo", type=Path, default=REPO)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    rel = args.release

    manifest = rel / "refs" / "manifest.tsv"
    releases = {}
    if manifest.exists():
        with open(manifest, newline="") as fh:       # descriptions hold quoted newlines
            for row in csv.DictReader(fh, delimiter="\t"):
                releases[row["uniprot_release"]] = sorted(
                    {*releases.get(row["uniprot_release"], []), row["date_fetched"]})
    explicit = Path("/opt/osmo_refdb/conda-explicit.txt")
    reuse_from = os.environ.get("REUSE_REFS_FROM", "")

    record = {
        "about": ("Build record of this release: which code, tools and inputs produced it, and "
                  "checksums of what came out. It does not change anything; it lets two builds "
                  "be compared (hpc/verify_reproducible.py) and a difference be traced to its "
                  "cause (a different commit, tool version, snapshot or machine)."),
        "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "release": rel.name,
        "code": {"osmo_refdb": git_info(args.repo), "osmotool": osmotool_info()},
        "environment": {
            "machine": platform.machine(), "cpu": cpu_model(), "host": platform.node(),
            "platform": platform.platform(), "python": sys.version.split()[0],
            "container": os.environ.get("APPTAINER_NAME") or os.environ.get("SINGULARITY_NAME")
                         or ("docker" if Path("/.dockerenv").exists() else "none"),
            "image_sha256": os.environ.get("SIF_SHA256", "not recorded"),
            "conda_explicit_list_sha256": sha256_of(explicit) if explicit.exists() else "not in image",
            "tools": tool_versions(),
        },
        "inputs": {
            "families_yaml_sha256": sha256_of(args.families) if args.families.exists() else "missing",
            "uniprot_releases_and_fetch_dates": releases or "no refs/manifest.tsv",
            "refs_reused_from": reuse_from or "none (fetched in this build)",
            "reused_alignments": os.environ.get("REUSE_ALIGNMENTS", "1" if reuse_from else "0"),
            "split_mode": os.environ.get("SPLIT_MODE", "random"),
            "refs": tree_hash(rel / "refs", "*.faa"),
            "alignments": tree_hash(rel / "alignments"),
        },
        "outputs": {
            "hmms_sha256_without_DATE_and_COM_lines": {
                p.name: hmm_sha256(p) for p in sorted((rel / "hmms").glob("*.hmm"))},
            "files": {n: sha256_of(rel / n) for n in (
                "osmo_refdb.dmnd", "osmo_refdb.diamond_cutoffs.tsv", "osmo_refdb.profile_cascade.tsv",
                "osmo_refdb.profile_report_levels.tsv", "osmo_refdb.profile_excluded_families.txt",
                "osmo_refdb.annotate_excluded_families.txt", "osmo_refdb.train_refs.faa")
                if (rel / n).exists()},
        },
    }
    catalogs = {}
    for cdir in sorted((rel / "catalogs").glob("*/")) if (rel / "catalogs").is_dir() else []:
        meta_path = cdir / "catalog.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        catalogs[cdir.name] = {
            "gene_file": meta.get("genes", ""), "gene_file_sha256": meta.get("genes_sha256", ""),
            "annotated_with_hmm_sha256": meta.get("hmm_sha256", ""),
            "index": meta.get("index", ""), "aligner": meta.get("aligner", ""),
            "n_genes_input": meta.get("n_genes_input"),
            "n_genes_with_hit": meta.get("n_genes_annotated_or_hit"),
            "annotation_sha256": sha256_of(cdir / "catalog.annotation.tsv")
                                 if (cdir / "catalog.annotation.tsv").exists() else "missing",
            "index_files": tree_hash(cdir, "catalog.fna.*")}
    record["outputs"]["catalogs"] = catalogs or "none built into this release"
    out = args.out or rel / "PROVENANCE.json"
    out.write_text(json.dumps(record, indent=2) + "\n")
    c = record["code"]["osmo_refdb"]
    print(f"Provenance written: {out}")
    print(f"  osmo_refdb {c.get('commit', '?')[:10]}{' (UNCOMMITTED CHANGES in the code)' if c.get('uncommitted_changes') else ''}, "
          f"osmotool {record['code']['osmotool'].get('version', '?')} "
          f"{'(shadowed by PYTHONPATH!)' if record['code']['osmotool'].get('shadowed_by_pythonpath') else ''}, "
          f"{record['environment']['machine']}")


if __name__ == "__main__":
    main()
