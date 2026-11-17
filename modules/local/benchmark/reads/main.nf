// DIAMOND (osmotool profile) and HMM (translate + hmmscan) on the simulated reads.
// mode = score: calibration run, every DIAMOND call kept with its HMM score
//      | on:    test run, with the cascade config decided on the calibration reads (cascade_in)
process BENCHMARK_READS {
    // calls: osmotool profile (diamond, hmmsearch), orfm, hmmscan
    container 'docker.io/barbarahelena/osmotool:0.8.0'
    label 'process_medium'
    tag "${mode}"
    input:
    path(reads)
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE): its DIAMOND files, cascade and hmms/, staged as db/
    val(mode)

    output:
    path('results'), emit: results
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('orfm'), eval("orfm --version 2>&1 | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_orfm

    script:
    """
    CASCADE=${mode} THREADS=${task.cpus} JOBS=${task.cpus} \\
        bash ${moduleDir}/templates/run_benchmark.sh ${reads} results db db/hmms/osmo_refdb.hmm
    """

    stub:
    """
    mkdir -p results/diamond results/hmm
    """
}
