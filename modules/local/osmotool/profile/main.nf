// osmotool profile: paired reads against an osmo_refdb release (DIAMOND, the DIAMOND+HMM cascade,
// copies-per-genome markers). Options go in ext.args (conf/modules.config): --keep_aln, --no_markers,
// --catalog, --no_catalog, filters.
process OSMOTOOL_PROFILE {
    tag "${meta.id}"
    label 'process_high'

    container 'docker.io/barbarahelena/osmotool:0.8.0'

    input:
    tuple val(meta), path(reads)   // [R1, R2]
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE): its DIAMOND files, cascade and hmms/, staged as db/

    output:
    tuple val(meta), path("${prefix}.gene_counts.tsv"),              emit: gene_counts
    tuple val(meta), path("${prefix}.blastx.tsv"),    optional: true, emit: alignments
    tuple val(meta), path("${prefix}.*"),                            emit: results
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('orfm'), eval("orfm --version 2>&1 | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_orfm

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    osmotool profile db -1 ${reads[0]} -2 ${reads[1]} \\
        --out_prefix ${prefix} --threads ${task.cpus} ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'gene\\traw_count\\n' > ${prefix}.gene_counts.tsv
    touch ${prefix}.blastx.tsv
    """
}
