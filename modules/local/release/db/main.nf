// The release as osmotool reads it (its DATABASE argument): one directory with the DIAMOND database
// and its files, hmms/, and the cascade config and report levels when there are any. Named after the
// release, since osmotool and the NaCl model record that name. Real files, not links: it is staged
// into every osmotool step. The rest of the release (refs, alignments, results) is published by the
// steps that make it (conf/modules.config).
process RELEASE_DB {
    tag "${name}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    val(name)
    path(hmms)
    path(release_files, stageAs: 'release_files/*')
    path(cascade, stageAs: 'cascade/*')               // [] for none: the calibration benchmark's database

    output:
    path("${name}"), emit: db

    script:
    """
    mkdir ${name}
    cp -L release_files/* ${name}/
    if [ -d cascade ]; then cp -L cascade/* ${name}/; fi
    cp -RL ${hmms} ${name}/hmms
    test -s ${name}/osmo_refdb.dmnd && test -s ${name}/hmms/osmo_refdb.hmm
    """

    stub:
    """
    mkdir -p ${name}/hmms
    touch ${name}/osmo_refdb.dmnd ${name}/hmms/osmo_refdb.hmm ${name}/hmms/osmo_refdb.hmm.h3m
    """
}
