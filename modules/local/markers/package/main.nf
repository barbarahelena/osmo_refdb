// The release's markers/: the universal single-copy marker genes' own DIAMOND database (built from
// markers.yaml by the MARKERS subworkflow, apart from the osmo families so they never compete with them),
// for osmotool profile's copies-per-genome normalisation (templates/package_markers.py).
// Published into the release (conf/modules.config).
process MARKERS_PACKAGE {
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(db_files, stageAs: 'marker_release/*')                  // the markers' DIAMOND_BUILD: dmnd, train refs, excluded lists
    path(fetch_manifest, stageAs: 'marker_release/refs/*')       // their REFS_COMBINE manifest.tsv: UniProt release, fetch date
    path(markers)                                                 // markers.yaml

    output:
    path('markers'), emit: markers
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python ${moduleDir}/templates/package_markers.py --marker-release marker_release --release . \\
        --markers ${markers} --out markers
    """

    stub:
    """
    mkdir markers && touch markers/osmo_markers.dmnd markers/osmo_markers.manifest.tsv
    """
}
