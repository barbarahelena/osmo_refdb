// Step 1 of a release's references: each family's positives and hard negatives from UniProt, and the HMM
// decoys' sequences (templates/fetch_refs.py; families.yaml says what to query). UniProt is live: the
// release and fetch dates are recorded in refs/manifest.tsv, and PROVENANCE keeps them.
// With --refs_store, kept there (storeDir, conf/modules.config) under <uniprot_release>/<key>, key a hash
// of the families file and fetch_refs.py: a later build in the same UniProt release reuses the fetch.
// So that it can be stored, refs/ is its one output and its tool versions are UNIPROT_RELEASE's.
process REFS_FETCH {
    tag "${families.baseName}"
    label 'process_medium'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:ddb8dc66d39a1fc1' :
        'community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:73160cfec41e18d3' }"

    input:
    tuple val(uniprot_release), val(key), path(families)   // UNIPROT_RELEASE; key: the stored fetch's

    output:
    path('refs'), emit: refs

    script:
    def T = "${moduleDir}/../templates"
    """
    python ${T}/fetch_refs.py --families ${families} --out refs --cache-dir cache
    # the fetch is of the release it is keyed on (UniProt can release in between)
    got=\$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "uniprot_release") c = i; next } c { print \$c }' refs/manifest.tsv | sort -u)
    if [ "\$got" != "${uniprot_release}" ]; then
        echo "ERROR: fetched UniProt release '\$got', expected ${uniprot_release}: rerun" >&2
        exit 1
    fi
    """

    stub:
    """
    mkdir refs
    printf '>a\\nMK\\n>b\\nMR\\n>c\\nMW\\n>d\\nMV\\n>e\\nMA\\n>f\\nMC\\n' > refs/ectA.positive.faa
    printf '>a\\nMK\\n>b\\nMR\\n' > refs/ectB.positive.faa
    printf '>n\\nMQ\\n' > refs/ectA.negative.faa
    printf 'family\\tuniprot_release\\tdate_fetched\\n' > refs/manifest.tsv
    printf 'ectA\\t${uniprot_release}\\ttoday\\n' >> refs/manifest.tsv
    """
}
