// The input set of each HMM decoy (families.yaml `hmm_decoys`): its fetched sequences, or a family's
// negative.train for a `from_negatives_of` decoy, minus every accession the benchmark could score
// (templates/build_hmm_decoys.py --select-only). The decoys are then clustered, aligned, trimmed and built
// by their own steps (subworkflows/local/release_build).
process HMM_DECOYSELECT {
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(refs)          // SPLIT_REFS: every family's sets, and the decoys' fetched sequences
    path(families)

    output:
    path('decoys/*.hmm_decoy.faa'), optional: true, emit: inputs
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    python ${moduleDir}/templates/build_hmm_decoys.py --select-only --refs ${refs} --families ${families} --out decoys
    """

    stub:
    """
    mkdir decoys
    printf '>a\\nMK\\n>b\\nMR\\n' > decoys/pepABC_decoy.hmm_decoy.faa
    """
}
