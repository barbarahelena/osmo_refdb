// Step 1 of a release's references, put together: every family's fetch (REFS_FETCHFAMILY) and the shared one
// (REFS_FETCHSHARED) in one refs/, their manifest rows in one refs/manifest.tsv. Checks that every family was
// fetched from the UniProt release the fetches are keyed on (UniProt can release between two of them).
process REFS_COMBINE {
    tag "${families.baseName}"
    label 'process_single'

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:bb86db01ccb6fea3' :
        'community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:8b81aff3d2d73d3c' }"

    input:
    val(uniprot_release)                                   // UNIPROT_RELEASE
    path(families)
    path(fetched, stageAs: 'fetched/*')                    // REFS_FETCHFAMILY: every family's directory
    path(shared)                                           // REFS_FETCHSHARED

    output:
    path('refs'), emit: refs

    script:
    """
    mkdir refs
    first=1
    for d in fetched/*/; do
        for f in "\$d"*; do
            [ "\$(basename "\$f")" = manifest.tsv ] && continue
            cp -L "\$f" refs/
        done
        if [ \$first = 1 ]; then cat "\$d"manifest.tsv; first=0; else tail -n +2 "\$d"manifest.tsv; fi >> refs/manifest.tsv
    done
    cp -L ${shared}/* refs/ 2>/dev/null || true
    # every family in families.yaml was fetched, from the release the fetches are keyed on
    python3 - <<'PY'
import csv, sys, yaml
want = [f["name"] for f in yaml.safe_load(open("${families}"))["families"]]
rows = list(csv.DictReader(open("refs/manifest.tsv"), delimiter="\\t"))
got = [r["family"] for r in rows]
missing = sorted(set(want) - set(got))
if missing or len(got) != len(set(got)):
    sys.exit(f"ERROR: fetched families do not match ${families}: missing {missing}, fetched {got}")
releases = sorted({r["uniprot_release"] for r in rows if r.get("uniprot_release")})
if releases != ["${uniprot_release}"]:
    sys.exit(f"ERROR: fetched UniProt release(s) {releases}, expected ${uniprot_release}: rerun")
order = {n: i for i, n in enumerate(want)}
rows.sort(key=lambda r: order[r["family"]])
with open("refs/manifest.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\\t")
    w.writeheader(); w.writerows(rows)
PY
    """

    stub:
    """
    mkdir refs
    for d in fetched/*/; do cp -L "\$d"*.faa refs/; done
    printf 'family\\tuniprot_release\\tdate_fetched\\n' > refs/manifest.tsv
    for d in fetched/*/; do tail -n +2 "\$d"manifest.tsv >> refs/manifest.tsv; done
    """
}
