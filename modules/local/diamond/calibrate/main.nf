// The DIAMOND database's release files besides the database itself: each family's bit-score cutoff
// from its calibration searches (templates/calibrate_diamond_cutoffs.py --hits), and the families
// kept out of profile's and annotate's output (templates/write_scope_manifest.py).
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
    path("${name}.{profile,annotate}_excluded_families.txt"),                             emit: excluded
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    mkdir -p hits
    python ${moduleDir}/templates/calibrate_diamond_cutoffs.py --refs ${refs} --release . \\
        --release-name ${name} --families ${families} --hits hits --out ${name}.diamond_cutoffs.tsv
    python ${moduleDir}/templates/write_scope_manifest.py --families ${families} \\
        --profile-out ${name}.profile_excluded_families.txt \\
        --annotate-out ${name}.annotate_excluded_families.txt
    """

    stub:
    """
    touch ${name}.diamond_cutoffs.tsv ${name}.profile_excluded_families.txt ${name}.annotate_excluded_families.txt
    """
}
