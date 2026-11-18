// Step 1 of a release's references: each family's positives and hard negatives from UniProt, and the HMM
// decoys' sequences (templates/fetch_refs.py; families.yaml says what to query). UniProt is live: the
// release and fetch dates are recorded in refs/manifest.tsv, and PROVENANCE keeps them.
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
    path(families)

    output:
    path('refs'),   emit: refs
    path('refs/*'), emit: files
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests

    script:
    def T = "${moduleDir}/../templates"
    """
    python ${T}/fetch_refs.py --families ${families} --out refs --cache-dir cache
    """

    stub:
    """
    mkdir refs
    printf '>a\\nMK\\n>b\\nMR\\n>c\\nMW\\n>d\\nMV\\n>e\\nMA\\n>f\\nMC\\n' > refs/ectA.positive.faa
    printf '>a\\nMK\\n>b\\nMR\\n' > refs/ectB.positive.faa
    printf '>n\\nMQ\\n' > refs/ectA.negative.faa
    printf 'family\\tuniprot_release\\tdate_fetched\\n' > refs/manifest.tsv
    """
}
