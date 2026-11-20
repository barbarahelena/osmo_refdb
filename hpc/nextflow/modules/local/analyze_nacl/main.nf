// The BacDive analysis on the new annotation, then the NaCl model.
process ANALYZE_NACL {
    // calls: python + numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, xgboost>=3.1, pyyaml
    // image: an image with those packages
    container 'todo.invalid/analyze_nacl'
    label 'process_high'
    input:
    val(release)
    val(annotated)

    output:
    val(release)

    script:
    """
    cd ${params.repo}
    PV=releases/${release}/phenotype_validation
    python pipeline/phenotype_validation/build_matrix.py --dir \$PV
    python pipeline/phenotype_validation/analyze.py --dir \$PV
    python pipeline/phenotype_validation/analyze_extended.py --dir \$PV
    if [ -d \$PV/strategy_control/annotate ]; then
        python pipeline/phenotype_validation/analyze_strategy_control.py --dir ${params.data_dir} --osmo-dir \$PV
    fi
    THREADS=${task.cpus} PYTHON=python bash run_pipeline.sh nacl_model ${release}
    """

    stub:
    """
    echo stub
    """
}
