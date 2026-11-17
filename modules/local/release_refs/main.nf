// The reference sequences of the new release: the source release's refs, minus its train/calib/test
// splits, minus study-tagged sequences, plus UniProt sequences for HMM decoys added to families.yaml since.
process RELEASE_REFS {
    // calls: python + requests, pyyaml; cd-hit (strip_ref_tag.py); network (UniProt)
    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:ddb8dc66d39a1fc1' :
        'community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:73160cfec41e18d3' }"
    label 'process_medium'
    input:
    path(families)
    path(source, stageAs: 'source')

    output:
    path('refs'), emit: refs
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests
    tuple val("${task.process}"), val('cd-hit'), eval("cd-hit -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_cd_hit

    script:
    """
    mkdir refs
    cp -R source/refs/. refs/
    rm -f refs/*.positive.train.faa refs/*.positive.calib.faa refs/*.positive.test.faa \\
          refs/*.negative.train.faa refs/*.negative.test.faa refs/*.decoy.faa refs/extra_positives_manifest.tsv
    python ${moduleDir}/templates/strip_ref_tag.py --refs refs --tag _study --threads ${task.cpus}
    python ${moduleDir}/templates/01_fetch_refs.py --families ${families} --out refs \\
        --only-hmm-decoys --cache-dir cache
    """

    stub:
    """
    mkdir refs
    """
}
