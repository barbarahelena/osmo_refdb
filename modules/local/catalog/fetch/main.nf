// The human-gut GMGC catalog, 2.99 GB out of a 19 GB Zenodo zip by range requests.
// storeDir: stored next to where the pipeline expects it; skipped if it is already there.
process CATALOG_FETCH {
    // calls: python (stdlib only), bash
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"
    label 'process_single'
    errorStrategy 'retry'   // downloads are resumable: a dropped connection is not a failure
    maxRetries 3
    storeDir "${file(params.catalog_genes).parent}"
    output:
    path('GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz'), emit: genes

    script:
    """
    python3 ${moduleDir}/templates/fetch_zenodo_member.py \\
        --member gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz \\
        --out GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
    test -s GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz

    """

    stub:
    """
    echo stub > GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
    """
}
