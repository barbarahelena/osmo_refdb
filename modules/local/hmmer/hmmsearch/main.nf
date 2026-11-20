// hmmsearch of one sequence set against one HMM: the per-sequence table (--tblout), <meta.id>.tblout.
// No --cut_ga: the scores are what calibration chooses the cutoff from.
process HMMER_HMMSEARCH {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--h7d74f8d_5' :
        'biocontainers/hmmer:3.4--h7d74f8d_5' }"

    input:
    tuple val(meta), path(hmm), path(seqs)

    output:
    tuple val(meta), path("${prefix}.tblout"), emit: tblout
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hmmsearch --cpu ${task.cpus} --tblout ${prefix}.tblout --noali ${args} ${hmm} ${seqs} > /dev/null
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tblout
    """
}
