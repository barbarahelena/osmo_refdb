// Annotate the gene catalog with the release's HMMs and package it into the release
// (annotation + the small hits-only index, catalogs/<name>/).
process ANNOTATE_CATALOG {
    // calls: osmotool build-catalog (hmmsearch, bwa-mem3), python for 14_package_catalog.py
    // image: osmotool (containers/README.md)
    container 'barbarahelena/osmotool:todo'
    label 'process_high'
    input:
    val(release)
    val(genes)

    output:
    val(release)

    script:
    """
    cd ${params.repo}
    python ${moduleDir}/templates/14_package_catalog.py --release releases/${release} \\
        --name ${params.catalog_name} --genes ${genes} --index hits --threads ${task.cpus} \\
        --source "${params.catalog_source}" --source-url "${params.catalog_source_url}" ${params.catalog_force ? '--force' : ''}
    test -s releases/${release}/catalogs/${params.catalog_name}/catalog.annotation.tsv
    """

    stub:
    """
    echo stub
    """
}
