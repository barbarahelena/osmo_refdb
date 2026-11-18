// Scores profile --catalog on the test reads, per reporting unit and family, in the three benchmark
// catalogs (A and B, indexed as the release's catalog is), against the same reads' translated calls (the test benchmark's
// DIAMOND run with the release's cascade): see templates/catalog_benchmark.py. The output is the
// release's results/catalog_benchmark/ (docs/catalog_benchmark.md).
process BENCHMARK_CATALOGSCORE {
    // calls: python + pandas, scikit-learn, matplotlib, pyyaml (compute_metrics.py)
    // Seqera container 'metrics' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scikit-learn:53d0d02acba1bb21' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_pyyaml_scikit-learn:55ad9ba88a0d577d' }"
    label 'process_single'
    input:
    path(reads)                                   // BENCHMARK_SIMULATE_TEST
    path(run, stageAs: 'run')                     // BENCHMARK_READS_TEST: the translated calls
    path(cascade, stageAs: 'cascade/*')           // cascade config and report levels
    path(families)
    path(runs, stageAs: 'catalog_runs/*')         // BENCHMARK_CATALOGMAP: one directory per catalog

    output:
    path('catalog_benchmark'), emit: metrics
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas

    script:
    """
    python ${moduleDir}/templates/catalog_benchmark.py --metrics-module ${moduleDir}/../metrics/templates/compute_metrics.py \\
        --reads ${reads} --diamond run/diamond --cascade cascade --families ${families} \\
        --runs catalog_runs --out catalog_benchmark
    """

    stub:
    """
    mkdir catalog_benchmark
    """
}
