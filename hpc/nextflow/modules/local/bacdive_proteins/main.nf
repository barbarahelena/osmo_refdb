// Release-independent, per genome: size/GC, Prodigal proteins (cached, gzipped) and the
// random-panel pool search, for the BacDive genomes and the salt-in control's.
// Task i of nchunks takes every nchunks-th genome.
process BACDIVE_PROTEINS {
    // calls: prodigal, hmmsearch, awk, find, xargs, sort, gzip
    // image: an image with: prodigal, hmmer, gawk, findutils, gzip (osmotool's image has prodigal and hmmer too)
    container 'todo.invalid/bacdive_proteins'
    label 'process_medium'
    tag "chunk ${chunk}/${params.nchunks}"
    input:
    val(chunk)
    val(data_dir)

    output:
    val(chunk)

    script:
    """
    cd ${params.repo}
    POOL=${data_dir}/random_panel/pool.hmm
    [ -s \$POOL ] || { echo "ERROR: \$POOL missing (snapshot, see hpc/README.md)" >&2; exit 1; }
    for PV in ${data_dir} ${data_dir}/strategy_control; do
        [ -d \$PV/genomes ] || continue
        CHUNK=${chunk}/${params.nchunks} bash pipeline/phenotype_validation/run_random_panel.sh \$PV \$POOL ${task.cpus}
    done
    """

    stub:
    """
    echo stub
    """
}
