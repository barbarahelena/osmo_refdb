// The NaCl-max model: trained on the BacDive strains' osmotool annotation (train.py), explained (SHAP,
// explain.py), and packaged as the release's nacl_model/ for osmotool (package_nacl_model.py).
// nacl_model_build/ (training table, CV, explain/) and nacl_model/ are published into the release
// (conf/modules.config).
process NACL_MODEL {
    tag "${release_name}"
    label 'process_high'

    // Seqera container 'nacl' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:92ae9dc030c0ee62' :
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:edc7132220048e3f' }"

    input:
    val(release_name)
    val(osmotool_version)                      // BACDIVE_ANNOTATE: the osmotool the training genomes were annotated with
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE): its DIAMOND files, cascade and hmms/, staged as db/
    path(matrix)                               // BACDIVE_MATRIX
    path(annotate, stageAs: 'annotate/*')      // BACDIVE_ANNOTATE: <bid>/ per strain
    path(strains, stageAs: 'strains/*')        // BACDIVE_PROTEINS: genome stats per strain
    path(family_ranking)                       // BACDIVE_ANALYZE

    output:
    path('nacl_model_build'), emit: build
    path('nacl_model'),       emit: model
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pandas'), eval("python3 -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('scikit-learn'), eval("python3 -c 'import sklearn; print(sklearn.__version__)'"), topic: versions, emit: versions_scikit_learn
    tuple val("${task.process}"), val('xgboost'), eval("python3 -c 'import xgboost; print(xgboost.__version__)'"), topic: versions, emit: versions_xgboost
    tuple val("${task.process}"), val('matplotlib'), eval("python3 -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib

    script:
    def T = "${moduleDir}/templates"
    """
    python ${T}/nacl_model/train.py --matrix ${matrix} --annotate annotate --strains strains \\
        --release db --release-name ${release_name} --osmotool-version ${osmotool_version} \\
        --out nacl_model_build --threads ${task.cpus}
    python ${T}/nacl_model/explain.py --build nacl_model_build --family-ranking ${family_ranking}
    python ${T}/package_nacl_model.py --build nacl_model_build --release db --release-name ${release_name} --out nacl_model
    """

    stub:
    """
    mkdir -p nacl_model_build/explain nacl_model
    """
}
