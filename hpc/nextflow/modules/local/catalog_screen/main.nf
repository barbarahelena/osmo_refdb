// THE GATE. Which families claim catalog genes that look like something else. The
// report is printed at the end of the run; a flag means "go and look".
process SCREEN {
    // calls: python + pyyaml (pipeline/15_catalog_screen.py)
    // image: an image with: python 3.11, pyyaml
    container 'todo.invalid/catalog_screen'
    label 'process_single'
    input:
    val(release)

    output:
    path('gate.txt'), emit: report

    script:
    """
    TASK_DIR=\$PWD
    cd ${params.repo}
    python pipeline/15_catalog_screen.py --release releases/${release} --catalog ${params.catalog_name} | tee \$TASK_DIR/screen.log
    {
        echo
        echo "================ THE GATE ================"
        grep '^Flagged:' \$TASK_DIR/screen.log
        echo "Read ${params.repo}/releases/${release}/catalogs/${params.catalog_name}/screen.tsv"
        echo "A flag means: look at those genes. A real confuser class means a decoy and a rebuild"
        echo "(new --release name). A false alarm means carry on: run the same command again with"
        echo "    -resume --gate_passed true"
        echo "=========================================="
    } > \$TASK_DIR/gate.txt
    """

    stub:
    """
    echo "stub gate" > gate.txt
    """
}
