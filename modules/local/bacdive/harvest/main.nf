// The full BacDive set (--bacdive_full): every strain record from the BacDive API (harvest_bacdive.sh
// sweep, about 1,820 calls, 1.5 h), the phenotype table of those with a NaCl maximum and a genome
// accession (build_table.py), and the within-domain salt-in control's strains (build_strategy_control.py).
// snapshot: harvested batches (bacdive_raw/, --bacdive_snapshot) that are not fetched again, so a
// shipped snapshot reproduces the published table (BacDive itself is live); [] harvests everything.
// The tables and batches are published to <repo>/<data_dir> (conf/modules.config).
process BACDIVE_HARVEST {
    label 'process_single'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // calls: curl, python (stdlib only), awk
    // Seqera container 'bacdive-fetch' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/curl_gawk_python:2d801c4a5165bf61' :
        'community.wave.seqera.io/library/curl_gawk_python:de4685943c4289cd' }"

    input:
    path(snapshot, stageAs: 'snapshot')

    output:
    path('bacdive/pilot_strains.tsv'),                    emit: strains          // the BacDive set
    path('bacdive/bacdive_classes.tsv'),                  emit: classes
    path('bacdive/strategy_control/pilot_strains.tsv'),   emit: control_strains  // the salt-in control
    path('bacdive/bacdive_raw'),                          emit: batches
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('curl'), eval("curl --version | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_curl

    script:
    def S = "${moduleDir}/templates/phenotype_validation"
    """
    mkdir -p bacdive/bacdive_raw
    if [ -d snapshot ]; then cp -RL snapshot/. bacdive/bacdive_raw/; fi
    bash ${S}/harvest_bacdive.sh bacdive sweep
    python3 ${S}/build_table.py --dir bacdive
    python3 ${S}/build_strategy_control.py --dir bacdive
    """

    stub:
    """
    mkdir -p bacdive/bacdive_raw bacdive/strategy_control
    printf 'bacdive_id\\tinsdc\\n1\\tGCA_000001\\n' > bacdive/pilot_strains.tsv
    printf 'bacdive_id\\ttax_class\\n' > bacdive/bacdive_classes.tsv
    printf 'bacdive_id\\tinsdc\\n' > bacdive/strategy_control/pilot_strains.tsv
    """
}
