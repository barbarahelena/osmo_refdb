// osmotool profile on one sample, with the catalog's full index if it exists. The reads are
// staged by Nextflow; the profile goes to ${params.outdir}/profile (publishDir, conf/modules.config).
process PROFILE_SAMPLE {
    // calls: osmotool profile (diamond, hmmsearch, orfm, bwa-mem3)
    // image: the osmotool image
    container 'todo.invalid/profile_sample'
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}*"), emit: profile

    script:
    """
    EXTRA=""
    [ -f ${params.full_index_dir}/catalog_index.json ] && EXTRA="--catalog_index ${params.full_index_dir}"
    osmotool profile ${params.repo}/releases/${params.release} -1 ${r1} -2 ${r2} \\
        --out_prefix ${meta.id} --threads ${task.cpus} \\
        --catalog ${params.catalog_name} \$EXTRA --keep_aln
    """

    stub:
    """
    echo stub > ${meta.id}.profile.tsv
    """
}
