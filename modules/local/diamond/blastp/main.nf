// diamond blastp: a protein set against a DIAMOND database, the hit table in <prefix>.tsv.
// ext.args (conf/modules.config): output format and thresholds.
process DIAMOND_BLASTP {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.8--he361c42_0' :
        'biocontainers/diamond:2.2.8--he361c42_0' }"

    input:
    tuple val(meta), path(query)
    path(db)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: hits
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    diamond blastp --query ${query} --db ${db} --out ${prefix}.tsv --threads ${task.cpus} ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
