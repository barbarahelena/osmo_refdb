// The human-gut GMGC catalog, 2.99 GB out of a 19 GB Zenodo zip by range requests.
// storeDir: stored next to where the pipeline expects it; skipped if it is already there.
process FETCH_GMGC {
    // calls: python (stdlib only), bash
    // image: an image with: python 3.11
    container 'todo.invalid/fetch_gmgc'
    label 'process_single'
    errorStrategy 'retry'   // downloads are resumable: a dropped connection is not a failure
    maxRetries 3
    storeDir "${file(params.catalog_genes).parent}"
    output:
    path('GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz')

    script:
    """
    bash ${params.repo}/hpc/00_fetch_gmgc.sh \$PWD
    test -s GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
    """

    stub:
    """
    echo stub > GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
    """
}
