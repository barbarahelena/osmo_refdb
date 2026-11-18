// Train / calibration / test split of the references (held out, no leakage).
process SPLIT_REFS {
    // calls: python + pyyaml
    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"
    label 'process_single'
    input:
    path(refs_in, stageAs: 'refs_in/*')   // REFS_BUILD: every refs/ file
    path(families)

    output:
    path('refs'), emit: refs
    path('refs/*'), emit: files   // the same, as files: the release's refs/ is these and the other steps' references
    path('refs/*.positive.train.faa'), emit: train   // per family: what its alignment and HMM are built from
    path('refs/*.faa'),                emit: faa     // every set, for the HMMs' calibration scoring
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    cp -RL refs_in refs
    python ${moduleDir}/templates/split_train_test.py --refs refs --families ${families} \\
        --split-mode ${params.split_mode}
    """

    stub:
    """
    cp -RL refs_in refs
    printf '>a\\nMK\\n>b\\nMR\\n' > refs/ectA.positive.train.faa
    printf '>a\\nMK\\n>c\\nMW\\n' > refs/ectB.positive.train.faa
    printf '>a\\nMK\\n' > refs/ectC.positive.train.faa
    printf '>p\\nMK\\n' > refs/ectA.positive.calib.faa
    printf '>n\\nMW\\n' > refs/ectA.negative.train.faa
    printf 'family\\tstatus\\n' > refs/negative_purity_manifest.tsv
    """
}
