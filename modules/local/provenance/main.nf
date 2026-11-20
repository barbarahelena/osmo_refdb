// The provenance receipt, written last, so it covers the catalog and the BacDive and NaCl results as
// well as the build. Every step runs in a container of its own, so
// the tool versions come from the steps themselves (the `versions` topic), not from this step's image.
// It reports none to the topic itself: it waits for the topic to close.
// PROVENANCE.json and software_versions.yml are published into the release (conf/modules.config).
process PROVENANCE {
    tag "${release_name}"
    label 'process_single'

    // calls: python + pyyaml, git (the pipeline's own checkout, projectDir)
    // Seqera container 'provenance' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/git_python_pyyaml:64763cd90b56dd64' :
        'community.wave.seqera.io/library/git_python_pyyaml:241a3ffe7fa635b5' }"

    input:
    val(release_name)
    path(refs, stageAs: 'refs/*')                       // the release's refs/
    path(alignments, stageAs: 'alignments/*')           // HMM_BUILD: every alignment, raw and trimmed
    path(hmms)                                          // HMM_CALIBRATE
    path(release_files, stageAs: 'files/*')             // DIAMOND_BUILD: the DIAMOND database and its files
    path(cascade, stageAs: 'files/*')                   // the calibrated cascade config and report levels
    tuple val(meta), path(catalog)                      // CATALOG_PACKAGE
    path(versions, stageAs: 'collected_versions.yml')   // software_versions.yml (the versions topic)
    path(families)
    path(pipeline)                                      // the pipeline's checkout: its git commit
    val(after)                                          // waits for the rest of this invocation's work

    output:
    path('PROVENANCE.json'),                              emit: provenance
    path('software_versions.yml'),                        emit: versions

    script:
    """
    cp -L ${versions} software_versions.yml
    SPLIT_MODE="${params.split_mode}" \\
    READ_LENGTHS="${params.read_lengths ?: ''}" BACKGROUND_R1="${params.background_r1 ?: ''}" BACKGROUND_R2="${params.background_r2 ?: ''}" \\
    DIAMOND_MIN_SUBJECT_COVER="${params.diamond_min_subject_cover}" \\
    python ${moduleDir}/templates/write_provenance.py --release-name ${release_name} \\
        --refs refs --alignments alignments --hmms ${hmms} --files files/* --catalogs ${catalog} \\
        --repo ${pipeline} --families ${families} --versions software_versions.yml --out PROVENANCE.json
    """

    stub:
    """
    echo '{}' > PROVENANCE.json
    cp -L ${versions} software_versions.yml
    """
}
