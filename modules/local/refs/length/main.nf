// Sequences whose length is far off their family's median (fusions, fragments) are flagged and dropped
// (templates/check_length_outliers.py; refs/length_outlier_manifest.tsv).
process REFS_LENGTH {
    tag "${families.baseName}"
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(refs_in, stageAs: 'refs_in')     // the previous step's refs/
    path(families)

    output:
    path('refs'),   emit: refs
    path('refs/*'), emit: files
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def T = "${moduleDir}/../templates"
    """
    cp -RL refs_in refs
    python ${T}/check_length_outliers.py --refs refs --families ${families}
    """

    stub:
    """
    cp -RL refs_in refs
    """
}
