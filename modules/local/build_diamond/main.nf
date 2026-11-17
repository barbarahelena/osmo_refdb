// DIAMOND database from the train references plus decoy and fused-ORF references, its per-family
// bit-score cutoffs, and the lists of families kept out of reported output.
process BUILD_DIAMOND {
    // calls: diamond, python + numpy, pyyaml
    // Seqera container 'diamond-py' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/diamond_numpy_python_pyyaml:5b3b1ed5ce857be9' :
        'community.wave.seqera.io/library/diamond_numpy_python_pyyaml:9cd53c8d486ed198' }"
    label 'process_medium'
    input:
    path(refs_in, stageAs: 'refs_in')
    path(families)

    output:
    path('refs'),          emit: refs
    path('release_files/*'), emit: release_files   // the DIAMOND database and its files: top level of the release
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    cp -R refs_in refs
    mkdir release_files
    T=${moduleDir}/templates
    python \$T/08a_build_decoy_refs.py --refs refs --families ${families}
    python \$T/08d_build_fusion_refs.py --refs refs --families ${families}
    THREADS=${task.cpus} bash \$T/08_build_diamond_db.sh refs release_files osmo_refdb
    python \$T/08b_calibrate_diamond_cutoffs.py --refs refs --release release_files \\
        --release-name osmo_refdb --families ${families} --threads ${task.cpus}
    python \$T/08c_write_scope_manifest.py --families ${families} \\
        --profile-out release_files/osmo_refdb.profile_excluded_families.txt \\
        --annotate-out release_files/osmo_refdb.annotate_excluded_families.txt
    """

    stub:
    """
    cp -R refs_in refs
    mkdir release_files && touch release_files/osmo_refdb.dmnd
    """
}
