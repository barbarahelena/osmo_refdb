// The random pool's candidate families: a seeded draw from Pfam-A, minus every Pfam accession in
// families.yaml and every family whose name or description matches an osmotic term
// (templates/build_random_pool.py candidates). ext.args (conf/modules.config): how many, the seed.
process POOL_CANDIDATES {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    path(pfam)
    path(families)

    output:
    path('candidates.hmm'), emit: hmm
    path('candidates.tsv'), emit: table
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def args = task.ext.args ?: ''
    """
    python ${moduleDir}/../templates/build_random_pool.py candidates --pfam ${pfam} --families ${families} --out-dir . ${args}
    """

    stub:
    """
    touch candidates.hmm
    printf 'pfam\\tname\\tdesc\\n' > candidates.tsv
    """
}
