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
    # every family in the families file gets its sets, so the stub follows whichever file the run uses
    python3 -c "import yaml; print('\\n'.join(f['name'] for f in yaml.safe_load(open('${families}'))['families']))" |
    while read -r fam; do
        printf '>a\\nMK\\n>b\\nMR\\n' > refs/\${fam}.positive.train.faa
        printf '>p\\nMK\\n' > refs/\${fam}.positive.calib.faa
        printf '>n\\nMW\\n' > refs/\${fam}.negative.train.faa
    done
    printf 'family\\tstatus\\n' > refs/negative_purity_manifest.tsv
    """
}
