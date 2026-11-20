// osmotool build-catalog: annotate a gene catalog with a release's HMMs and index it (bwa-mem3).
// ext.args (conf/modules.config): --index neighbours|hits|all|none; ext.prefix: the output directory's name.
process OSMOTOOL_BUILDCATALOG {
    tag "${meta.id}"
    label 'process_high'

    container 'docker.io/barbarahelena/osmotool:0.9.0'

    input:
    tuple val(meta), path(genes)
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE): its DIAMOND files, cascade and hmms/, staged as db/

    output:
    tuple val(meta), path("${prefix}"), emit: catalog
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('bwa-mem3'), eval("bwa-mem3 version 2>&1 | head -1"), topic: versions, emit: versions_bwa_mem3

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    osmotool build-catalog db --genes ${genes} -o ${prefix} --threads ${task.cpus} ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    echo '{}' > ${prefix}/catalog.json
    touch ${prefix}/catalog.annotation.tsv
    """
}
