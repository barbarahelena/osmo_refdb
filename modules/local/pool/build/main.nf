// The random pool: the candidates found in at least two screening proteomes, drawn stratified by
// prevalence so each tenth-of-prevalence bin has ext.args' --per-family times as many pool families as
// the osmo panel has there (in the pilot's feature matrix) (templates/build_random_pool.py pool).
// pool.hmm and pool_families.tsv are what BACDIVE_PROTEINS searches and BACDIVE_ANALYZE draws from;
// published to <repo>/<data_dir>/random_panel (conf/modules.config).
process POOL_BUILD {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    path(pfam)
    path(candidates)            // POOL_CANDIDATES: candidates.tsv
    path(screen)                // the candidates' hmmsearch --cut_ga against the screening proteomes
    path(osmo_matrix)           // the screening genomes' osmo feature matrix
    val(n_genomes)

    output:
    path('pool.hmm'),          emit: hmm
    path('pool_families.tsv'), emit: families
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def args = task.ext.args ?: ''
    """
    python ${moduleDir}/../templates/build_random_pool.py pool --pfam ${pfam} --out-dir . \\
        --screen ${screen} --osmo-matrix ${osmo_matrix} --n-genomes ${n_genomes} ${args}
    """

    stub:
    """
    printf 'HMMER3/f\\nNAME  PF00001\\n//\\n' > pool.hmm
    printf 'pfam\\tname\\tdesc\\tscreen_genomes_hit\\n' > pool_families.tsv
    """
}
