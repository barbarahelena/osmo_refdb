// Annotate the gene catalog with the release's HMMs and package it into the release
// (annotation + the small hits-only index, catalogs/<name>/).
process ANNOTATE_CATALOG {
    // calls: osmotool build-catalog (hmmsearch, bwa-mem3), python for pipeline/14_package_catalog.py
    // image: the osmotool image (a release with catalog mode)
    container 'todo.invalid/annotate_catalog'
    label 'process_high'
    input:
    val(release)
    val(genes)

    output:
    val(release)

    script:
    """
    cd ${params.repo}
    CATALOG_NAME=${params.catalog_name} CATALOG_GENES=${genes} CATALOG_INDEX=hits \\
    CATALOG_SOURCE="${params.catalog_source}" CATALOG_SOURCE_URL="${params.catalog_source_url}" \\
    CATALOG_FORCE=${params.catalog_force ? 1 : 0} THREADS=${task.cpus} \\
        bash run_pipeline.sh catalog ${release}
    test -s releases/${release}/catalogs/${params.catalog_name}/catalog.annotation.tsv
    """

    stub:
    """
    echo stub
    """
}
