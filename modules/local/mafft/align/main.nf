// MAFFT multiple sequence alignment of one sequence set. ext.args (conf/modules.config): --auto.
// One thread by default (process_single): MAFFT with several threads is not reproducible run to run.
process MAFFT_ALIGN {
    tag "${meta.id}"
    label 'process_single'

    // Seqera container 'align' (environment.yml next to this file): mafft 7.526 is not on biocontainers
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/trimal_mafft_python_pyyaml:0aafcb62c9cac1a6' :
        'community.wave.seqera.io/library/trimal_mafft_python_pyyaml:d3ab7e2b4cc1b639' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.aln.fasta"), emit: alignment
    tuple val("${task.process}"), val('mafft'), eval("mafft --version 2>&1 | grep -oE 'v[0-9.]+' | head -1"), topic: versions, emit: versions_mafft

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mafft ${args} --thread ${task.cpus} ${fasta} > ${prefix}.aln.fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln.fasta
    """
}
