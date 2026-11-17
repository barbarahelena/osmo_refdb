// Train / calibration / test split of the references (held out, no leakage).
process SPLIT_REFS {
    // calls: python + pyyaml
    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"
    label 'process_single'
    input:
    path(refs_in, stageAs: 'refs_in')
    path(families)

    output:
    path('refs'), emit: refs
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    cp -R refs_in refs
    python ${moduleDir}/templates/03_split_train_test.py --refs refs --families ${families} \\
        --split-mode ${params.split_mode}
    """

    stub:
    """
    cp -R refs_in refs
    """
}
