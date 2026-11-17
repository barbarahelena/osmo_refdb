// trimAl: trims an alignment. ext.args (conf/modules.config): the gap threshold, -gt, per family
// from families.yaml (meta.trim_gt).
process TRIMAL {
    tag "${meta.id}"
    label 'process_single'

    // Seqera container 'align' (environment.yml next to this file): the same build as MAFFT_ALIGN's
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/trimal_mafft_python_pyyaml:0aafcb62c9cac1a6' :
        'community.wave.seqera.io/library/trimal_mafft_python_pyyaml:d3ab7e2b4cc1b639' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${prefix}.trimmed.fasta"), emit: trimmed
    tuple val("${task.process}"), val('trimal'), eval("trimal --version 2>&1 | grep -oE 'v[0-9][0-9.a-z]*' | head -1 | sed 's/^v//'"), topic: versions, emit: versions_trimal

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    trimal -in ${alignment} -out ${prefix}.trimmed.fasta ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.trimmed.fasta
    """
}
