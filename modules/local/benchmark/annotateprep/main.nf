// The protein set for the annotate benchmark: held-out positives, hard negatives and background,
// each on its own one-gene contig with a Prodigal-style header (benchmark/templates/benchmark_annotate.py, shared with the other annotate benchmark step; --prepare).
process BENCHMARK_ANNOTATEPREP {
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(refs, stageAs: 'refs/*')
    path(families)

    output:
    tuple path('annotate/dummy_assembly.fna'), path('annotate/eval_proteins.faa'), emit: eval_set
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    python ${moduleDir}/../templates/benchmark_annotate.py --prepare --refs refs --families ${families} --out annotate
    """

    stub:
    """
    mkdir annotate && touch annotate/dummy_assembly.fna annotate/eval_proteins.faa
    """
}
