// osmotool index-catalog: the full, competitive index of a gene catalog. Depends only on the genes
// and the aligner, not on a release. ext.args: --aligner bwa (classic bwa; bwa-mem3 needs 100+ GB
// for a whole gut catalog); ext.prefix: the index directory's name.
// The tool versions are written into the index directory (versions.yml), not emitted: under storeDir
// the index is kept outside the work directory and the step is skipped once it exists.
// NEEDS an osmotool with catalog mode: not in 0.8.0 (it is on osmotool's catalog-profile branch).
process OSMOTOOL_INDEXCATALOG {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/barbarahelena/osmotool:0.8.0'

    input:
    tuple val(meta), path(genes)

    output:
    tuple val(meta), path("${prefix}"), emit: index

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    osmotool index-catalog ${args} --genes ${genes} -o ${prefix}
    test -f ${prefix}/catalog_index.json

    cat <<-END_VERSIONS > ${prefix}/versions.yml
    "${task.process}":
        osmotool: \$(osmotool --version | sed 's/^osmotool //')
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    echo '{}' > ${prefix}/catalog_index.json
    printf '"%s":\\n    osmotool: stub\\n' "${task.process}" > ${prefix}/versions.yml
    """
}
