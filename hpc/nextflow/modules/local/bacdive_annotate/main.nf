// osmotool annotate with the release's HMMs, per genome, from the cached proteins.
process BACDIVE_ANNOTATE {
    // calls: osmotool annotate (hmmsearch, diamond, prodigal)
    // image: the osmotool image
    container 'todo.invalid/bacdive_annotate'
    label 'process_low'
    tag "chunk ${chunk}/${params.nchunks}"
    input:
    val(chunk)
    val(release)

    output:
    val(chunk)

    script:
    """
    cd ${params.repo}
    PV=releases/${release}/phenotype_validation
    for D in \$PV \$PV/strategy_control; do
        [ -d \$D/genomes ] || continue
        CHUNK=${chunk}/${params.nchunks} bash pipeline/phenotype_validation/run_annotate_cached.sh \$D releases/${release} ${task.cpus}
    done
    """

    stub:
    """
    echo stub
    """
}
