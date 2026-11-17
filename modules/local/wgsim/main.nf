// wgsim: paired reads from a genome. The read names carry the true coordinates.
// ext.args (conf/modules.config): seed, number of pairs, read lengths.
process WGSIM {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wgsim:1.0--h577a1d6_10' :
        'biocontainers/wgsim:1.0--h577a1d6_10' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_R{1,2}.fastq"), emit: reads
    tuple val("${task.process}"), val('wgsim'), eval("wgsim 2>&1 | sed -n 's/^Version: *//p'"), topic: versions, emit: versions_wgsim

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    wgsim ${args} ${fasta} ${prefix}_R1.fastq ${prefix}_R2.fastq > /dev/null
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.fastq ${prefix}_R2.fastq
    """
}
