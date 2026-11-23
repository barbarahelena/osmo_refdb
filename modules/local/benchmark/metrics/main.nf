// Precision / recall / F1 for DIAMOND, DIAMOND+cascade and HMM.
// It imports ../../templates/report_rule.py (the report-level rule, shared with the annotate benchmark).
// mode = calib: scores the calibration reads and WRITES the cascade config and report levels (cascade/),
//               with metrics/decision_sensitivity.tsv: the decisions on point estimates alone vs the support rule
//      | test:  scores the test reads with those (cascade_in)
process BENCHMARK_METRICS {
    // calls: python + numpy, pandas, matplotlib, pyyaml
    // Seqera container 'metrics' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scikit-learn:53d0d02acba1bb21' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scikit-learn:55ad9ba88a0d577d' }"
    label 'process_single'
    tag "${mode}"
    input:
    path(reads)
    path(run, stageAs: 'run')
    path(cascade_in, stageAs: 'cascade_in/*')   // mode test: the calibration's cascade files
    path(families)
    val(mode)

    output:
    path('metrics'),                     emit: metrics
    path('cascade/*'), optional: true,   emit: cascade   // mode calib: cascade config and report levels
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('scikit-learn'), eval("python3 -c 'import sklearn; print(sklearn.__version__)'"), topic: versions, emit: versions_scikit_learn
    tuple val("${task.process}"), val('matplotlib'), eval("python3 -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib

    script:
    def args = mode == 'calib'
        ? '--cascade-out cascade/osmo_refdb.profile_cascade.tsv --report-levels-out cascade/osmo_refdb.profile_report_levels.tsv'
        : '--cascade-config cascade_in/osmo_refdb.profile_cascade.tsv --report-levels cascade_in/osmo_refdb.profile_report_levels.tsv'
    def extra = task.ext.args ?: ''   // e.g. --min-eval-positives (conf/modules.config)
    """
    mkdir cascade
    python ${moduleDir}/templates/compute_metrics.py --reads ${reads} --diamond run/diamond --hmm run/hmm \\
        --out metrics --families ${families} ${args} ${extra}
    """

    stub:
    """
    mkdir metrics cascade
    if [ "${mode}" = calib ]; then
        touch cascade/osmo_refdb.profile_cascade.tsv cascade/osmo_refdb.profile_report_levels.tsv metrics/decision_sensitivity.tsv
    fi
    """
}
