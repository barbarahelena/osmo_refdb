// hmmbuild: a profile HMM from a (trimmed) protein alignment, named after meta.id.
process HMMER_HMMBUILD {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--h7d74f8d_5' :
        'biocontainers/hmmer:3.4--h7d74f8d_5' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${prefix}.hmm"), emit: hmm
    tuple val("${task.process}"), val('hmmer'), eval("hmmbuild -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hmmbuild --amino -n ${meta.id} --cpu ${task.cpus} ${args} ${prefix}.hmm ${alignment} > /dev/null
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'HMMER3/f\\nNAME  ${meta.id}\\nLENG  10\\n//\\n' > ${prefix}.hmm
    """
}
