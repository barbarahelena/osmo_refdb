// BacDive strains and genomes into data_dir. Resumable; batches already in
// <data_dir>/bacdive_raw are skipped, so a shipped snapshot there reproduces the table.
process FETCH_BACDIVE {
    // calls: curl, python (stdlib only), awk, gzip, xargs, sed, grep
    // image: an image with: curl, python 3.11, gawk, gzip, findutils, sed, grep
    container 'todo.invalid/fetch_bacdive'
    label 'process_high'
    errorStrategy 'retry'   // downloads are resumable: a dropped connection is not a failure
    maxRetries 3
    output:
    val(params.data_dir)

    script:
    """
    cd ${params.repo}
    PAR=${task.cpus} bash hpc/00_fetch_bacdive.sh ${params.data_dir}
    test -d ${params.data_dir}/genomes
    """

    stub:
    """
    echo stub
    """
}
