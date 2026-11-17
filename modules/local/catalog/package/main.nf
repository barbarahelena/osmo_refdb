// Packages the catalog built by OSMOTOOL_BUILDCATALOG as the release's catalogs/<name>/: provenance in
// catalog.json, the hmmsearch evidence gzipped, SHA256SUMS (templates/package_catalog.py).
// Published into releases/<release>/catalogs/ (conf/modules.config).
process CATALOG_PACKAGE {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    val(release_name)
    tuple val(meta), path(catalog, stageAs: 'built/*')
    path(genes)

    output:
    tuple val(meta), path("${meta.id}"), emit: catalog
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python ${moduleDir}/templates/package_catalog.py --release-name ${release_name} \\
        --name ${meta.id} --genes ${genes} --built ${catalog} --out ${meta.id} \\
        --source "${params.catalog_source}" --source-url "${params.catalog_source_url}"
    test -s ${meta.id}/catalog.annotation.tsv
    """

    stub:
    """
    mkdir ${meta.id} && touch ${meta.id}/catalog.annotation.tsv && echo '{}' > ${meta.id}/catalog.json
    """
}
