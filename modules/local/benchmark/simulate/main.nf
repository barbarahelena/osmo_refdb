// Simulated read pairs from one split of the references, with per-read truth from CDS overlap.
// split = calib (positive.calib + negative.train: used to decide cascade and report levels)
//       | test  (held-out positives, hard negatives, random background: used to score)
process BENCHMARK_SIMULATE {
    // calls: wgsim, python + numpy, pyyaml
    // Seqera container 'simulate' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/insilicoseq_wgsim_python_pyyaml:a4bccdbc78a6f5da' :
        'community.wave.seqera.io/library/insilicoseq_wgsim_python_pyyaml:e5b214e4c49c8d0d' }"
    label 'process_low'
    tag "${split}"
    input:
    path(refs, stageAs: 'refs/*')
    path(families)
    val(split)

    output:
    path('reads'), emit: reads
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml
    tuple val("${task.process}"), val('iss'), eval("iss --version 2>&1 | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_iss
    tuple val("${task.process}"), val('wgsim'), eval("wgsim 2>&1 | sed -n 's/^Version: *//p'"), topic: versions, emit: versions_wgsim

    script:
    """
    python ${moduleDir}/templates/09_simulate_reads.py --refs refs --families ${families} \\
        --reads-per-sequence 20 --split ${split} --out reads
    python ${moduleDir}/templates/09b_compute_read_truth.py --reads reads
    """

    stub:
    """
    mkdir reads
    """
}
