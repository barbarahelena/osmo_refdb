// For the families that opt in (families.yaml refseq_gene_symbols), more positives from NCBI RefSeq
// (templates/add_refseq_positives.py; refs/refseq_positives_manifest.tsv).
process REFS_REFSEQ {
    tag "${families.baseName}"
    label 'process_single'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:ddb8dc66d39a1fc1' :
        'community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:73160cfec41e18d3' }"

    input:
    path(refs_in, stageAs: 'refs_in')     // the previous step's refs/
    path(families)

    output:
    path('refs'),   emit: refs
    path('refs/*'), emit: files
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests

    script:
    def T = "${moduleDir}/../templates"
    """
    cp -RL refs_in refs
    python ${T}/add_refseq_positives.py --refs refs --families ${families}
    """

    stub:
    """
    cp -RL refs_in refs
    """
}
