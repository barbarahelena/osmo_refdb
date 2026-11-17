// Release-independent, per genome: size/GC, Prodigal proteins (gzipped, for BACDIVE_ANNOTATE) and the
// random-panel pool search, for the BacDive genomes and the salt-in control's. Task i of nchunks takes
// every nchunks-th genome; each strain's results are one directory, <bid>/.
// Prodigal and hmmsearch are osmotool's own: annotate --proteins reads these proteins in its place.
process BACDIVE_PROTEINS {
    tag "chunk ${chunk}/${params.nchunks}"
    label 'process_medium'

    // calls: prodigal, hmmsearch, awk, find, xargs, sort, gzip
    container 'docker.io/barbarahelena/osmotool:0.8.0'

    input:
    val(chunk)
    path(genomes)                                       // BACDIVE_FETCH
    path(control_genomes, stageAs: 'control_genomes')   // BACDIVE_FETCH: the salt-in control's
    path(pool)                                          // the random-panel pool, pool.hmm

    output:
    tuple val(chunk), path('strains'), path('control_strains'), emit: chunk   // for BACDIVE_ANNOTATE's chunk
    path('strains/*'),         optional: true, emit: strains                    // per strain, for the analysis
    path('control_strains/*'), optional: true, emit: control_strains
    tuple val("${task.process}"), val('prodigal'), eval("prodigal -v 2>&1 | grep -oE 'V[0-9.]+' | tr -d V"), topic: versions, emit: versions_prodigal
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    def run = "${moduleDir}/templates/phenotype_validation/run_random_panel.sh"
    """
    export CHUNK=${chunk}/${params.nchunks}
    bash ${run} ${genomes} strains ${pool} ${task.cpus}
    bash ${run} control_genomes control_strains ${pool} ${task.cpus}
    """

    stub:
    """
    mkdir -p strains/${chunk}00 control_strains
    touch strains/${chunk}00/${chunk}00.faa.gz
    """
}
