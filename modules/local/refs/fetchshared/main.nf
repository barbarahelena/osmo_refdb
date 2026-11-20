// Step 1 of a release's references, what needs every family once (templates/fetch_refs.py --shared-only): the
// dedicated fusion-pair candidates (from the families' fetched positives and Pfam domains) and the HMM decoys'
// sequences. With --refs_store kept there (storeDir) under <uniprot_release>/<families>.shared_<key>, key a hash
// of families.yaml's hmm_decoys, its fusion_partner families' entries and fetch_refs.py.
process REFS_FETCHSHARED {
    tag "${families.baseName}"
    label 'process_low'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:bb86db01ccb6fea3' :
        'community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:8b81aff3d2d73d3c' }"

    input:
    tuple val(uniprot_release), val(key), path(families)   // UNIPROT_RELEASE; key: the stored fetch's
    path(fetched, stageAs: 'fetched/*')                    // REFS_FETCHFAMILY: every family's directory

    output:
    path('shared'), emit: refs

    script:
    def T = "${moduleDir}/../templates"
    """
    mkdir -p evidence shared
    for d in fetched/*/; do ln -s "\$(cd "\$d" && pwd)"/* evidence/; done
    python ${T}/fetch_refs.py --families ${families} --out shared --cache-dir cache \\
        --shared-only --evidence-dir evidence
    """

    stub:
    """
    mkdir shared
    """
}
