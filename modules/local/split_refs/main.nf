// Train / calibration / test split of the references (held out, no leakage).
process SPLIT_REFS {
    // calls: python + pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_single'
    input:
    path(refs_in, stageAs: 'refs_in')
    path(families)

    output:
    path('refs'), emit: refs

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
