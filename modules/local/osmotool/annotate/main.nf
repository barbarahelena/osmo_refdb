// osmotool annotate: an assembly (or pre-called proteins) against an osmo_refdb release.
// Options go in ext.args (conf/modules.config): --method, --keep_proteins, --keep_aln, filters.
process OSMOTOOL_ANNOTATE {
    tag "${meta.id}"
    label 'process_low'

    container 'docker.io/barbarahelena/osmotool:0.12.5'

    input:
    tuple val(meta), path(fasta), path(proteins)   // proteins: [] to let osmotool run Prodigal
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE): its DIAMOND files, cascade and hmms/, staged as db/

    output:
    tuple val(meta), path("${prefix}.gene_counts.tsv"),                     emit: gene_counts
    tuple val(meta), path("${prefix}.systems.tsv"),          optional: true, emit: systems
    tuple val(meta), path("${prefix}.gene_coordinates.tsv"), optional: true, emit: coordinates
    tuple val(meta), path("${prefix}.prodigal.faa"),         optional: true, emit: proteins
    tuple val(meta), path("${prefix}.*"),                                   emit: results
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('prodigal'), eval("prodigal -v 2>&1 | grep -oE 'V[0-9.]+' | tr -d V"), topic: versions, emit: versions_prodigal

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def prot = proteins ? "--proteins ${proteins}" : ''
    """
    osmotool annotate db ${fasta} ${prot} \\
        --out_prefix ${prefix} --threads ${task.cpus} ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'gene\\traw_count\\n' > ${prefix}.gene_counts.tsv
    printf 'system\\tpresent\\tsubunit_count\\trescued_hits\\n' > ${prefix}.systems.tsv
    printf 'protein_id\\tgene\\n' > ${prefix}.gene_coordinates.tsv
    echo ">${meta.id}_1 # 1 # 300 # 1 # stub" > ${prefix}.prodigal.faa
    """
}
