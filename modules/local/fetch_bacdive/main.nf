// BacDive strains and genomes into data_dir. Resumable; batches already in
// <data_dir>/bacdive_raw are skipped, so a shipped snapshot there reproduces the table.
process FETCH_BACDIVE {
    // calls: curl, python (stdlib only), awk, gzip, xargs, sed, grep
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_high'
    errorStrategy 'retry'   // downloads are resumable: a dropped connection is not a failure
    maxRetries 3
    output:
    val(params.data_dir)

    script:
    """
    SCRIPTS=${moduleDir}/templates/phenotype_validation
    PV=${params.repo}/${params.data_dir}
    mkdir -p \$PV
    cd ${params.repo}
    # 1. harvest (batches already in bacdive_raw are skipped: a shipped snapshot there reproduces the table)
    bash \$SCRIPTS/harvest_bacdive.sh \$PV sweep
    # 2. phenotype table, 3. genomes from NCBI
    python3 \$SCRIPTS/build_table.py --dir \$PV
    PAR=${task.cpus} bash \$SCRIPTS/fetch_genomes.sh \$PV
    # 4. the within-domain salt-in control: strains and genomes
    python3 \$SCRIPTS/build_strategy_control.py --dir \$PV
    PAR=4 bash \$SCRIPTS/fetch_genomes.sh \$PV/strategy_control
    test -d \$PV/genomes
    """

    stub:
    """
    echo stub
    """
}
