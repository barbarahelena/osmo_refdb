// The BacDive analysis on the new annotation, then the NaCl model.
process ANALYZE_NACL {
    // calls: python + numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, xgboost>=3.1, pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
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
    python ${moduleDir}/templates/phenotype_validation/build_matrix.py --dir \$PV
    python ${moduleDir}/templates/phenotype_validation/analyze.py --dir \$PV
    python ${moduleDir}/templates/phenotype_validation/analyze_extended.py --dir \$PV
    if [ -d \$PV/strategy_control/annotate ]; then
        python ${moduleDir}/templates/phenotype_validation/analyze_strategy_control.py --dir ${params.data_dir} --osmo-dir \$PV
    fi
    T=${moduleDir}/templates
    python \$T/nacl_model/train.py --pv-dir \$PV --release releases/${release} \\
        --out releases/${release}/nacl_model_build --threads ${task.cpus}
    python \$T/nacl_model/explain.py --build releases/${release}/nacl_model_build --pv-dir \$PV
    python \$T/13_package_nacl_model.py --build releases/${release}/nacl_model_build --release releases/${release}
    """

    stub:
    """
    echo stub
    """
}
