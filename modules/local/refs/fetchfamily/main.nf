// A release's references (REFS_BUILD), per family: its positives and hard negatives from UniProt
// (templates/fetch_refs.py --skip-shared; families.yaml says what to query). UniProt is live: the release and
// fetch dates go into the family's manifest.tsv row, and PROVENANCE keeps them.
// One task per family, at most a few at a time (maxForks: UniProt answers 429 to more). With --refs_store,
// each family's fetch is kept there (storeDir, conf/modules.config) under <uniprot_release>/<name>_<key>, key a
// hash of that family's families.yaml entry and fetch_refs.py: a later build in the same UniProt release
// fetches only the families that changed. So that it can be stored, the family's directory is its one
// output, and its tool versions are UNIPROT_RELEASE's.
process REFS_FETCHFAMILY {
    tag "${name}"
    label 'process_low'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3
    maxForks 4

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:bb86db01ccb6fea3' :
        'community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:8b81aff3d2d73d3c' }"

    input:
    tuple val(uniprot_release), val(name), val(key), val(entry)   // UNIPROT_RELEASE; the family; key: its stored fetch's; entry: its families.yaml entry, as YAML

    output:
    path("${name}"), emit: refs

    script:
    def T = "${moduleDir}/../templates"
    """
    cat > family.yaml <<'FAMILY_YAML'
${entry}
FAMILY_YAML
    python ${T}/fetch_refs.py --families family.yaml --out ${name} --cache-dir cache --skip-shared
    # fetched from the release it is stored under: fetch_refs.py records "unknown" when its own release
    # lookup fails, and a stored fetch is never redone. Failing here retries it and keeps it out of the store.
    python3 - <<'PY'
import csv, sys
got = sorted({r["uniprot_release"] for r in csv.DictReader(open("${name}/manifest.tsv"), delimiter="\\t")})
if got != ["${uniprot_release}"]:
    sys.exit(f"ERROR: ${name} fetched from UniProt release(s) {got}, expected ${uniprot_release}")
PY
    """

    stub:
    """
    mkdir ${name}
    printf '>a\\nMK\\n>b\\nMR\\n>c\\nMW\\n>d\\nMV\\n>e\\nMA\\n>f\\nMC\\n' > ${name}/${name}.positive.faa
    printf '>n\\nMQ\\n' > ${name}/${name}.negative.faa
    printf 'family\\tuniprot_release\\tdate_fetched\\n${name}\\t${uniprot_release}\\ttoday\\n' > ${name}/manifest.tsv
    """
}
