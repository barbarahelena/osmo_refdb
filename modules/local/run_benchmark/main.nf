// DIAMOND (osmotool profile) and HMM (translate + hmmscan) on the simulated reads.
// mode = score: calibration run, every DIAMOND call kept with its HMM score
//      | on:    test run, with the cascade config decided on the calibration reads (cascade_in)
process RUN_BENCHMARK {
    // calls: osmotool profile (diamond, hmmsearch), orfm, hmmscan
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_medium'
    tag "${mode}"
    input:
    path(reads)
    path(hmms)
    path(release_files)
    path(cascade_in, stageAs: 'cascade_in')
    val(mode)

    output:
    path('results'), emit: results

    script:
    """
    mkdir rel
    for f in release_files/* cascade_in/*; do [ -e \$f ] && ln -s \$PWD/\$f rel/\$(basename \$f); done
    ln -s \$PWD/hmms rel/hmms
    CASCADE=${mode} THREADS=${task.cpus} JOBS=${task.cpus} \\
        bash ${moduleDir}/templates/10_run_benchmark.sh reads results rel hmms/osmo_refdb.hmm
    """

    stub:
    """
    mkdir -p results/diamond results/hmm
    """
}
