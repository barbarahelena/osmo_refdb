// OrfM: open reading frames of reads (R1 and R2 together), for the HMM pass on reads.
process ORFM {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orfm:2.1.1--h54198d6_1' :
        'biocontainers/orfm:2.1.1--h54198d6_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.orfs.faa"), emit: orfs
    tuple val("${task.process}"), val('orfm'), eval("orfm --version 2>&1 | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_orfm

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    zcat -f ${reads} > merged.fastq
    orfm ${args} merged.fastq > ${prefix}.orfs.faa
    rm merged.fastq
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.orfs.faa
    """
}
