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
    refs/manifest.tsv, the pipeline parameters that shape the benchmark (from
    the environment: SPLIT_MODE, READ_LENGTHS, BACKGROUND_R1/R2,
    DIAMOND_MIN_SUBJECT_COVER), the parameters of the report-level and
    cascade decision, set per mode (DECISION_RULE_ENV: PROFILE_* and
    ANNOTATE_*, and the shared DECISION_BOOTSTRAP_SEED), and a
    checksum over refs/ and alignments/;
  * outputs: checksums of every HMM (without the DATE and COM lines HMMER
    writes), the DIAMOND database, cutoffs, cascade config and report levels;
  * licence: the code is MIT, the release is CC-BY-4.0 with attribution to
    the sources it is built from (UniProt; GMGC and the prevalence-filtered
    catalogs when a catalog is in it).

This is a record, not a guarantee: it says what ran so a rebuild can be
compared and a difference traced.

Usage:
  python bin/write_provenance.py --release releases/v10 [--families families.yaml]
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

REPO = Path(__file__).resolve().parents[1]

# the report-level and cascade rule's parameters, per mode (params.<name lowercased>)
DECISION_RULE_ENV = (
    "PROFILE_MIN_REPORT_F1", "PROFILE_MIN_EVAL_POSITIVES", "PROFILE_MIN_CASCADE_GAIN",
    "PROFILE_DECISION_SUPPORT", "PROFILE_DECISION_BOOTSTRAP",
    "ANNOTATE_MIN_REPORT_F1", "ANNOTATE_MIN_EVAL_POSITIVES", "ANNOTATE_DECISION_SUPPORT",
    "ANNOTATE_DECISION_BOOTSTRAP", "DECISION_BOOTSTRAP_SEED",
)

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


def process_versions(path: Path) -> dict:
    """The versions.yml every pipeline step writes (nf-core style), collected by the workflow:
    {process: {tool: version}}. Each step runs in a container of its own, so this, not a probe
    of this step's own image, is what the steps ran with."""
    import yaml
    out: dict = {}
    for doc in yaml.safe_load_all(path.read_text()):
        for process, tools in (doc or {}).items():
            out[process.split(":")[-1]] = {k: str(v) for k, v in (tools or {}).items()}
    return dict(sorted(out.items()))


def merged_tools(per_process: dict) -> dict:
    """{tool: version}, or {tool: [versions]} where steps disagree."""
    seen: dict = {}
    for tools in per_process.values():
        for k, v in tools.items():
            seen.setdefault(k, set()).add(v)
    return {k: (sorted(v)[0] if len(v) == 1 else sorted(v)) for k, v in sorted(seen.items())}


def cpu_model() -> str:
    try:
        for line in open("/proc/cpuinfo"):
            if line.lower().startswith(("model name", "hardware", "cpu model")):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def licence_record(has_catalogs: bool) -> dict:
    """The code's licence, the release's licence and whom the release must credit. They differ: the
    code is MIT, the release (databases, HMMs, reference sequences, gene catalogs) is CC-BY-4.0 and
    carries its sources' attribution. The catalog's sources only apply when a catalog is in the release."""
    sources = [{
        "source": "UniProt",
        "covers": "the reference sequences, and so the HMMs, the DIAMOND database and the benchmark sets",
        "cite": "The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. "
                "Nucleic Acids Res 53(D1):D609-D617 (2025)",
        "doi": "10.1093/nar/gkae1010",
        "licence": "CC-BY-4.0"}]
    if has_catalogs:
        sources += [{
            "source": "GMGC (Global Microbial Gene Catalog)",
            "covers": "catalogs/: the genes, their annotation and index",
            "cite": "Coelho LP et al. Towards the biogeography of prokaryotic genes. Nature 601:252-256 (2022)",
            "doi": "10.1038/s41586-021-04233-4"}, {
            "source": "Prevalence-filtered GMGC catalogs (Cayman)",
            "covers": "catalogs/: the filtered gene set the catalog is taken from",
            "cite": "Ducarmon QR et al. Cayman enables large-scale analysis of gut microbiome "
                    "carbohydrate-active enzyme repertoires. Nat Microbiol 11:1739-1753 (2026)",
            "doi": "10.1038/s41564-026-02318-2",
            "data": "https://zenodo.org/records/10473258",
            "licence": "CC-BY-4.0"}]
    return {
        "code": {"licence": "MIT", "file": "LICENSE", "repository": "https://github.com/barbarahelena/osmo_refdb"},
        "release": {
            "licence": "CC-BY-4.0",
            "url": "https://creativecommons.org/licenses/by/4.0/",
            "attribution": "Reuse of the release needs attribution to the sources below.",
            "sources": sources},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", type=Path, default=None,
                    help="a release directory; or give its parts with the options below")
    ap.add_argument("--release-name", default=None, help="default: the --release directory's name")
    ap.add_argument("--refs", type=Path, default=None, help="default: <release>/refs")
    ap.add_argument("--alignments", type=Path, default=None, help="default: <release>/alignments")
    ap.add_argument("--hmms", type=Path, default=None, help="default: <release>/hmms")
    ap.add_argument("--files", type=Path, nargs="*", default=None,
                    help="the release's top-level files (DIAMOND database, cutoffs, cascade, ...); "
                         "default: those in <release>")
    ap.add_argument("--catalogs", type=Path, nargs="*", default=None,
                    help="catalog directories; default: <release>/catalogs/*/")
    ap.add_argument("--families", type=Path, default=REPO / "families.yaml")
    ap.add_argument("--repo", type=Path, default=REPO)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--versions", type=Path, default=None,
                    help="the pipeline's collected versions.yml: tool versions per step, "
                         "instead of probing the tools in this environment")
    args = ap.parse_args()
    rel = args.release
    if rel is None and not (args.release_name and args.refs and args.hmms):
        ap.error("give --release, or --release-name, --refs and --hmms (and the other parts)")
    refs_dir = args.refs or rel / "refs"
    aln_dir = args.alignments or (rel / "alignments" if rel else None)
    hmm_dir = args.hmms or rel / "hmms"
    top = {f.name: f for f in args.files} if args.files is not None else \
        {f.name: f for f in rel.iterdir() if f.is_file()} if rel else {}
    catalog_dirs = args.catalogs if args.catalogs is not None else \
        sorted((rel / "catalogs").glob("*/")) if rel and (rel / "catalogs").is_dir() else []
    release_name = args.release_name or rel.name

    manifest = refs_dir / "manifest.tsv"
    releases = {}
    if manifest.exists():
        with open(manifest, newline="") as fh:       # descriptions hold quoted newlines
            for row in csv.DictReader(fh, delimiter="\t"):
                releases[row["uniprot_release"]] = sorted(
                    {*releases.get(row["uniprot_release"], []), row["date_fetched"]})
    per_process = process_versions(args.versions) if args.versions else {}
    osmotool = osmotool_info()
    if per_process and "error" in osmotool:     # not importable here: each step has its own image
        versions = merged_tools(per_process).get("osmotool")
        osmotool = {"version": versions or "not used", "source": "versions.yml of the steps",
                    # what `osmotool --version` printed in each step's container (a list above: they differ)
                    "per_step": {proc: tools["osmotool"] for proc, tools in sorted(per_process.items())
                                 if "osmotool" in tools}}
    explicit = Path("/opt/osmo_refdb/conda-explicit.txt")

    record = {
        "about": ("Build record of this release: which code, tools and inputs produced it, and "
                  "checksums of what came out. It does not change anything; it lets two builds "
                  "be compared and a difference be traced to its "
                  "cause (a different commit, tool version, snapshot or machine)."),
        "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "release": release_name,
        "code": {"osmo_refdb": git_info(args.repo), "osmotool": osmotool},
        "environment": {
            "machine": platform.machine(), "cpu": cpu_model(), "host": platform.node(),
            "platform": platform.platform(), "python": sys.version.split()[0],
            "container": os.environ.get("APPTAINER_NAME") or os.environ.get("SINGULARITY_NAME")
                         or ("docker" if Path("/.dockerenv").exists() else "none"),
            "image_sha256": os.environ.get("SIF_SHA256", "not recorded"),
            "conda_explicit_list_sha256": sha256_of(explicit) if explicit.exists() else "not in image",
            "tools": merged_tools(per_process) if per_process else tool_versions(),
            "tools_per_step": per_process or "not recorded (no --versions)",
        },
        "inputs": {
            "families_yaml_sha256": sha256_of(args.families) if args.families.exists() else "missing",
            "uniprot_releases_and_fetch_dates": releases or "no refs/manifest.tsv",
            "split_mode": os.environ.get("SPLIT_MODE", "random"),
            "benchmark_read_lengths": os.environ.get("READ_LENGTHS") or "simulator default",
            "benchmark_background_reads": [x for x in (os.environ.get("BACKGROUND_R1"), os.environ.get("BACKGROUND_R2")) if x] or "none",
            "benchmark_diamond_min_subject_cover": os.environ.get("DIAMOND_MIN_SUBJECT_COVER", "0.0"),
            # how the report levels and the cascade were decided, per mode (report_rule.py)
            "decision_rule": {k.lower(): os.environ.get(k, "not recorded") for k in DECISION_RULE_ENV},
            "refs": tree_hash(refs_dir, "*.faa"),
            "alignments": tree_hash(aln_dir) if aln_dir else {"n_files": 0},
        },
        "outputs": {
            "hmms_sha256_without_DATE_and_COM_lines": {
                p.name: hmm_sha256(p) for p in sorted(hmm_dir.glob("*.hmm"))},
            "files": {n: sha256_of(top[n]) for n in (
                "osmo_refdb.dmnd", "osmo_refdb.diamond_cutoffs.tsv", "osmo_refdb.profile_cascade.tsv",
                "osmo_refdb.profile_report_levels.tsv", "osmo_refdb.annotate_report_levels.tsv",
                "osmo_refdb.profile_excluded_families.txt",
                "osmo_refdb.annotate_excluded_families.txt", "osmo_refdb.train_refs.faa")
                if n in top},
        },
    }
    catalogs = {}
    for cdir in catalog_dirs:
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
    record["licence"] = licence_record(has_catalogs=bool(catalogs))
    out = args.out or (rel / "PROVENANCE.json" if rel else Path("PROVENANCE.json"))
    out.write_text(json.dumps(record, indent=2) + "\n")
    c = record["code"]["osmo_refdb"]
    print(f"Provenance written: {out}")
    print(f"  osmo_refdb {c.get('commit', '?')[:10]}{' (UNCOMMITTED CHANGES in the code)' if c.get('uncommitted_changes') else ''}, "
          f"osmotool {record['code']['osmotool'].get('version', '?')} "
          f"{'(shadowed by PYTHONPATH!)' if record['code']['osmotool'].get('shadowed_by_pythonpath') else ''}, "
          f"{record['environment']['machine']}")


if __name__ == "__main__":
    main()
