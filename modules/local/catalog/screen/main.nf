// The catalog screen: which families claim catalog genes that look like something else. The report is printed
// at the end of the run; a flag means "go and look". screen.tsv is published next to the catalog in
// the release (conf/modules.config).
process CATALOG_SCREEN {
    tag "${meta.id}"
    label 'process_single'

    // calls: python + pyyaml (catalog_screen.py)
    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(hmms)                          // HMM_CALIBRATE: the calibration scores, hmms/scores/
    tuple val(meta), path(catalog)      // CATALOG_PACKAGE: catalog.annotation.tsv
    path(families)                      // the HMM decoys

    output:
    tuple val(meta), path('screen.tsv'), emit: screen
    path('screen_report.txt'),           emit: report
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    python ${moduleDir}/templates/catalog_screen.py --scores-dir ${hmms}/scores \\
        --annotation ${catalog}/catalog.annotation.tsv --families ${families} --out screen.tsv | tee screen.log
    {
        echo
        echo "============ THE CATALOG SCREEN ============"
        grep '^Flagged:' screen.log || echo "Flagged: none"
        echo "In ${params.repo}/releases/${params.release}/catalogs/${meta.id}/screen.tsv; read it before"
        echo "publishing the release. A flag means: look at those genes. A real confuser class means a"
        echo "decoy and a rebuild (new --release name); a false alarm, nothing to do."
        echo "============================================"
    } > screen_report.txt
    """

    stub:
    """
    touch screen.tsv
    echo "stub screen" > screen_report.txt
    """
}
