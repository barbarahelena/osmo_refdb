// The benchmark's random bacterial background: 5000 UniProt proteins with none of the panel's Pfam
// domains, minus those >= 70% identical to a panel positive (DIAMOND) (templates/build_random_background.py;
// refs/random_background.faa). Benchmark only: never in a database or a calibration.
// NEEDS the 'refs' image rebuilt with diamond (environment.yml lists it).
process REFS_BACKGROUND {
    tag "${families.baseName}"
    label 'process_medium'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:ddb8dc66d39a1fc1' :
        'community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:73160cfec41e18d3' }"

    input:
    path(refs_in, stageAs: 'refs_in')     // the previous step's refs/
    path(families)

    output:
    path('refs'),   emit: refs
    path('refs/*'), emit: files
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests

    script:
    def T = "${moduleDir}/../templates"
    """
    cp -RL refs_in refs
    python ${T}/build_random_background.py --refs refs --families ${families} --threads ${task.cpus}
    """

    stub:
    """
    cp -RL refs_in refs
    """
}
