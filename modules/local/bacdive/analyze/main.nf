// The BacDive phenotype validation: does the osmo panel predict salt tolerance (analyze.py), better than
// random panels of the same size, with family ranking and cross-validation (analyze_extended.py), and
// the within-domain salt-in control (analyze_strategy_control.py). Each test's printed report is kept
// as <script>.txt. Published to releases/<release>/phenotype_validation (conf/modules.config).
process BACDIVE_ANALYZE {
    label 'process_medium'

    // Seqera container 'nacl' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:92ae9dc030c0ee62' :
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:edc7132220048e3f' }"

    input:
    path(matrix)                                                // BACDIVE_MATRIX
    path(pool_families)                                         // the random-panel pool, pool_families.tsv
    path(strains, stageAs: 'strains/*')                         // BACDIVE_PROTEINS: <bid>/ per strain
    path(classes)                                               // BACDIVE_FETCH: bacdive_classes.tsv
    path(control_pilot, stageAs: 'control_pilot_strains.tsv')   // BACDIVE_FETCH: the salt-in control
    path(control_strains, stageAs: 'control_strains/*')         // BACDIVE_PROTEINS: the control's
    path(control_annotate, stageAs: 'control_annotate/*')       // BACDIVE_ANNOTATE: the control's

    output:
    path('family_ranking.tsv'),                                 emit: family_ranking
    path('{null_panels,cv_results}.tsv'),                       emit: tables
    path('strategy_control/*.tsv'), optional: true,             emit: strategy_control
    path('*.txt'),                                              emit: reports
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('scipy'), eval("python3 -c 'import scipy; print(scipy.__version__)'"), topic: versions, emit: versions_scipy
    tuple val("${task.process}"), val('scikit-learn'), eval("python3 -c 'import sklearn; print(sklearn.__version__)'"), topic: versions, emit: versions_scikit_learn
    tuple val("${task.process}"), val('statsmodels'), eval("python3 -c 'import statsmodels; print(statsmodels.__version__)'"), topic: versions, emit: versions_statsmodels

    script:
    def T = "${moduleDir}/templates"
    // the control is only analysed when it has strains with annotation
    def control = control_annotate ? """
    python ${T}/analyze_strategy_control.py --matrix ${matrix} --classes ${classes} --strains strains \\
        --control-pilot ${control_pilot} --control-strains control_strains --control-annotate control_annotate \\
        --out-dir strategy_control | tee analyze_strategy_control.txt""" : ''
    """
    mkdir -p strains control_strains control_annotate
    python ${T}/analyze.py --matrix ${matrix} | tee analyze.txt
    python ${T}/analyze_extended.py --matrix ${matrix} --pool-families ${pool_families} --strains strains \\
        --out-dir . | tee analyze_extended.txt
    ${control}
    """

    stub:
    """
    touch family_ranking.tsv null_panels.tsv cv_results.tsv analyze.txt
    """
}
