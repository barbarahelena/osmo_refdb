// Each family's DIAMOND bit-score cutoff, from its calibration searches
// (templates/calibrate_diamond_cutoffs.py --hits): <name>.diamond_cutoffs.tsv.
process DIAMOND_CALIBRATE {
    tag "${name}"
    label 'process_single'

    // Seqera container 'diamond-py' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/diamond_numpy_python_pyyaml:5b3b1ed5ce857be9' :
        'community.wave.seqera.io/library/diamond_numpy_python_pyyaml:9cd53c8d486ed198' }"

    input:
    val(name)
    path(refs)                          // SPLIT_REFS: the calibration sets the searches were run on
    path(hits, stageAs: 'hits/*')       // DIAMOND_BLASTP: <family>.{positive,negative}.tsv
    path(families)

    output:
    path("${name}.diamond_cutoffs.tsv"),                                                  emit: cutoffs
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    mkdir -p hits
    python ${moduleDir}/templates/calibrate_diamond_cutoffs.py --refs ${refs} --release . \\
        --release-name ${name} --families ${families} --hits hits --out ${name}.diamond_cutoffs.tsv
    """

    stub:
    """
    touch ${name}.diamond_cutoffs.tsv
    """
}
