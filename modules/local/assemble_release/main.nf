// Put the pieces together as releases/<release> in the clone, where everything after the build
// (catalog, BacDive, validation, provenance) reads it. Refuses to overwrite an existing release.
// The simulated reads are not kept; the metrics are.
process ASSEMBLE_RELEASE {
    // calls: cp
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_single'
    input:
    path(refs)
    path(alignments)
    path(hmms)
    path(release_files)
    path(qc)
    path(cascade)
    path(calib_metrics, stageAs: 'calib_metrics')
    path(metrics)
    path(annotate)

    output:
    val(params.release)

    script:
    """
    D=${params.repo}/releases/${params.release}
    test ! -e \$D || { echo "ERROR: \$D exists; use a new release name or remove it" >&2; exit 1; }
    mkdir -p \$D/results
    cp -RL refs alignments hmms \$D/
    cp -RL release_files/. cascade/. \$D/
    cp -L qc_scorecard.tsv \$D/
    cp -RL calib_metrics \$D/results/calib_metrics
    cp -RL metrics annotate \$D/results/
    test -s \$D/hmms/osmo_refdb.hmm && test -s \$D/osmo_refdb.dmnd
    """

    stub:
    """
    echo stub
    """
}
