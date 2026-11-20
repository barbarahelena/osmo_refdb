// Scores osmotool annotate's calls on the benchmark protein set, per family and paralog group
// (benchmark/templates/benchmark_annotate.py, shared with the other annotate benchmark step; --calls-from). Writes annotate/ as the release keeps it in results/.
process BENCHMARK_ANNOTATESCORE {
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(refs, stageAs: 'refs/*')
    path(families)
    path(annotations, stageAs: 'calls/*')   // annotate_<method>.* from OSMOTOOL_ANNOTATE

    output:
    path('annotate'),     emit: annotate
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    python ${moduleDir}/../templates/benchmark_annotate.py --refs refs --families ${families} --out annotate --calls-from calls
    cp -L calls/* annotate/
    """

    stub:
    """
    mkdir annotate
    """
}
