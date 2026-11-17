// BacDive strains and genomes: harvest the BacDive API, the phenotype table, the genomes from NCBI,
// and the within-domain salt-in control (its strains, and genomes for those the main set lacks).
// snapshot: a bacdive_raw/ directory of harvested batches (--bacdive_snapshot); batches in it are not
// fetched again, so a shipped snapshot reproduces the published table (BacDive itself is live). [] to
// harvest everything. The tables and the batches are published to <repo>/<data_dir> (conf/modules.config).
process BACDIVE_FETCH {
    label 'process_high'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // calls: curl, python (stdlib only), awk, gzip, xargs, sed, grep
    // Seqera container 'bacdive-fetch' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/curl_gawk_python:2d801c4a5165bf61' :
        'community.wave.seqera.io/library/curl_gawk_python:de4685943c4289cd' }"

    input:
    path(snapshot, stageAs: 'snapshot')

    output:
    path('bacdive/genomes'),                              emit: genomes
    path('bacdive/pilot_strains.tsv'),                    emit: pilot
    path('bacdive/bacdive_classes.tsv'),                  emit: classes
    path('bacdive/strategy_control/genomes'),             emit: control_genomes
    path('bacdive/strategy_control/pilot_strains.tsv'),   emit: control_pilot
    path('bacdive/bacdive_raw'),                          emit: batches
    path('bacdive/genome_fetch.log'), optional: true,     emit: log
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('curl'), eval("curl --version | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_curl

    script:
    def S = "${moduleDir}/templates/phenotype_validation"
    """
    mkdir -p bacdive/bacdive_raw
    if [ -d snapshot ]; then cp -RL snapshot/. bacdive/bacdive_raw/; fi
    # 1. harvest (batches already there are skipped), 2. phenotype table, 3. genomes from NCBI
    bash ${S}/harvest_bacdive.sh bacdive sweep
    python3 ${S}/build_table.py --dir bacdive
    PAR=${task.cpus} bash ${S}/fetch_genomes.sh bacdive
    # 4. the within-domain salt-in control: strains, and genomes (links to the main set's where it has them)
    python3 ${S}/build_strategy_control.py --dir bacdive
    PAR=4 bash ${S}/fetch_genomes.sh bacdive/strategy_control
    test -d bacdive/genomes
    """

    stub:
    """
    mkdir -p bacdive/genomes bacdive/strategy_control/genomes bacdive/bacdive_raw
    echo '>c' > bacdive/genomes/1_GCA_000001.fna
    printf 'bacdive_id\\n1\\n' > bacdive/pilot_strains.tsv
    printf 'bacdive_id\\ttax_class\\n' > bacdive/bacdive_classes.tsv
    printf 'bacdive_id\\n' > bacdive/strategy_control/pilot_strains.tsv
    """
}
