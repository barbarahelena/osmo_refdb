// Precision / recall / F1 for DIAMOND, DIAMOND+cascade and HMM.
// mode = calib: scores the calibration reads and WRITES the cascade config and report levels (cascade/)
//      | test:  scores the test reads with those (cascade_in)
process COMPUTE_METRICS {
    // calls: python + numpy, pandas, matplotlib, pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_single'
    tag "${mode}"
    input:
    path(reads)
    path(run, stageAs: 'run')
    path(cascade_in, stageAs: 'cascade_in')
    path(families)
    val(mode)

    output:
    path('metrics'),           emit: metrics
    path('cascade'), optional: true, emit: cascade

    script:
    def args = mode == 'calib'
        ? '--cascade-out cascade/osmo_refdb.profile_cascade.tsv --report-levels-out cascade/osmo_refdb.profile_report_levels.tsv'
        : '--cascade-config cascade_in/osmo_refdb.profile_cascade.tsv --report-levels cascade_in/osmo_refdb.profile_report_levels.tsv'
    """
    mkdir results cascade
    ln -s \$PWD/reads results/reads; ln -s \$PWD/run/diamond results/diamond; ln -s \$PWD/run/hmm results/hmm
    python ${moduleDir}/templates/11_compute_metrics.py --results results --out metrics \\
        --families ${families} ${args}
    """

    stub:
    """
    mkdir metrics cascade
    """
}
