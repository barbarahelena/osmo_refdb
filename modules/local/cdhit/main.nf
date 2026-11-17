// CD-HIT: clusters a protein set, keeping one representative per cluster (<prefix>.faa, and the
// clusters in <prefix>.faa.clstr). ext.args (conf/modules.config): identity threshold and the rest.
process CDHIT {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5ca1c30_13' :
        'biocontainers/cd-hit:4.8.1--h5ca1c30_13' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.faa"),       emit: clustered
    tuple val(meta), path("${prefix}.faa.clstr"), emit: clusters
    tuple val("${task.process}"), val('cd-hit'), eval("cd-hit -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_cd_hit

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.clustered"
    """
    cd-hit -i ${fasta} -o ${prefix}.faa -T ${task.cpus} ${args} > /dev/null
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.clustered"
    """
    printf '>a\\nMK\\n>b\\nMR\\n' > ${prefix}.faa
    touch ${prefix}.faa.clstr
    """
}
