// Hard negatives that are true positives under another gene symbol (>= 70% identity to the family's
// positives, DIAMOND) are flagged and dropped (templates/check_negative_purity.py;
// refs/negative_purity_manifest.tsv).
process REFS_PURITY {
    tag "${families.baseName}"
    label 'process_medium'

    // Seqera container 'diamond-py' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/diamond_numpy_python_pyyaml:5b3b1ed5ce857be9' :
        'community.wave.seqera.io/library/diamond_numpy_python_pyyaml:9cd53c8d486ed198' }"

    input:
    path(refs_in, stageAs: 'refs_in')     // the previous step's refs/
    path(families)

    output:
    path('refs'),   emit: refs
    path('refs/*'), emit: files
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond

    script:
    def T = "${moduleDir}/../templates"
    """
    cp -RL refs_in refs
    python ${T}/check_negative_purity.py --refs refs --families ${families} --threads ${task.cpus}
    """

    stub:
    """
    cp -RL refs_in refs
    """
}
