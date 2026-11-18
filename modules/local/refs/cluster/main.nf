// One family's positives or negatives clustered at 90% identity (CD-HIT): redundancy removed, and no
// near-duplicates across the train/test split. Fewer than 5 sequences are kept as they are. Writes both
// <family>.<label>.clustered.faa and <family>.<label>.faa (the set the later steps use).
process REFS_CLUSTER {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cd-hit:4.8.1--h5ca1c30_13' :
        'biocontainers/cd-hit:4.8.1--h5ca1c30_13' }"

    input:
    tuple val(meta), path(fasta, stageAs: 'input.faa')   // meta.id: <family>.<positive|negative>

    output:
    tuple val(meta), path("${meta.id}.faa"), path("${meta.id}.clustered.faa"), emit: clustered
    tuple val(meta), path("${meta.id}.clustered.faa.clstr"), optional: true, emit: clusters
    tuple val("${task.process}"), val('cd-hit'), eval("cd-hit -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_cd_hit

    script:
    def args = task.ext.args ?: ''
    """
    n=\$(grep -c '^>' input.faa)
    if [ "\$n" -lt 5 ]; then
        echo "[${meta.id}] only \$n sequences: kept as they are"
        cp input.faa ${meta.id}.clustered.faa
    else
        cd-hit -i input.faa -o ${meta.id}.clustered.faa -T ${task.cpus} ${args} > /dev/null
        echo "[${meta.id}] \$n -> \$(grep -c '^>' ${meta.id}.clustered.faa) sequences"
    fi
    cp ${meta.id}.clustered.faa ${meta.id}.faa
    """

    stub:
    """
    cp input.faa ${meta.id}.clustered.faa
    cp input.faa ${meta.id}.faa
    """
}
