// osmotool annotate with the release's HMMs, per genome, from the proteins BACDIVE_PROTEINS cached (no
// Prodigal rerun), for the BacDive genomes and the salt-in control's. Task i takes the same genomes as
// BACDIVE_PROTEINS task i. Each strain's annotation is one directory, <bid>/.
process BACDIVE_ANNOTATE {
    tag "chunk ${chunk}/${params.nchunks}"
    label 'process_low'

    // calls: osmotool annotate (hmmsearch, diamond)
    container 'docker.io/barbarahelena/osmotool:0.8.0'

    input:
    tuple val(chunk), path(strains), path(control_strains)   // BACDIVE_PROTEINS, the same chunk
    path(genomes)                                              // BACDIVE_FETCH
    path(control_genomes, stageAs: 'control_genomes')
    path(db)                                                   // RELEASE_DB

    output:
    path('annotate/*'),         optional: true, emit: annotate
    path('control_annotate/*'), optional: true, emit: control_annotate
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    def run = "${moduleDir}/templates/phenotype_validation/run_annotate_cached.sh"
    """
    export CHUNK=${chunk}/${params.nchunks}
    bash ${run} ${genomes} ${strains} annotate ${db} ${task.cpus}
    bash ${run} control_genomes ${control_strains} control_annotate ${db} ${task.cpus}
    """

    stub:
    """
    mkdir -p annotate/${chunk}00 control_annotate
    touch annotate/${chunk}00/${chunk}00.gene_counts.tsv
    """
}
