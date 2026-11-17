// The strain x family matrix the BacDive analysis and the NaCl model start from: each strain's
// phenotype (pilot_strains.tsv) and its osmo panel gene counts (osmotool annotate), for strains with
// an assembly of at least 1 Mb.
process BACDIVE_MATRIX {
    label 'process_single'

    // Seqera container 'nacl' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:92ae9dc030c0ee62' :
        'community.wave.seqera.io/library/matplotlib_numpy_pandas_python_pruned:edc7132220048e3f' }"

    input:
    path(genomes)                                 // BACDIVE_FETCH
    path(pilot)                                   // BACDIVE_FETCH: pilot_strains.tsv
    path(annotate, stageAs: 'annotate/*')         // BACDIVE_ANNOTATE: <bid>/ per strain

    output:
    path('feature_matrix.tsv'), emit: matrix
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python ${moduleDir}/templates/build_matrix.py --genomes ${genomes} --pilot ${pilot} \\
        --annotate annotate --out feature_matrix.tsv
    """

    stub:
    """
    touch feature_matrix.tsv
    """
}
