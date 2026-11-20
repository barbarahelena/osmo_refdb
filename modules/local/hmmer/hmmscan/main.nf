// hmmscan of sequences against a pressed HMM database: hmm_dir holds <hmm_name> and its .h3[fimp]
// files. No --cut_ga by default: short read ORFs never reach cutoffs calibrated on full-length proteins.
process HMMER_HMMSCAN {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--h7d74f8d_5' :
        'biocontainers/hmmer:3.4--h7d74f8d_5' }"

    input:
    tuple val(meta), path(seqs)
    path(hmm_dir)
    val(hmm_name)

    output:
    tuple val(meta), path("${prefix}.hmmscan.tblout"), emit: tblout
    tuple val("${task.process}"), val('hmmer'), eval("hmmscan -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hmmscan --cpu ${task.cpus} --tblout ${prefix}.hmmscan.tblout --noali ${args} ${hmm_dir}/${hmm_name} ${seqs} > /dev/null
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hmmscan.tblout
    """
}
