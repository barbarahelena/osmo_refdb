// Scores osmotool annotate's calls on one split's protein set, per family and paralog group
// (benchmark/templates/benchmark_annotate.py, shared with the other annotate benchmark step; --calls-from).
// split = calib: also decides what annotate reports, levels/osmo_refdb.annotate_report_levels.tsv
//       | test:  also average precision per family, from the alignments osmotool kept (--keep_aln)
// Writes annotate/ as the release keeps it in results/ (the raw alignments stay in the work directory).
process BENCHMARK_ANNOTATESCORE {
    tag "${split}"
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(refs, stageAs: 'refs/*')
    path(families)
    path(annotations, stageAs: 'calls/*')   // annotate_<method>.* from OSMOTOOL_ANNOTATE
    val(split)                              // calib | test

    output:
    path('annotate'),                    emit: annotate
    path('levels/*'), optional: true,    emit: levels     // split calib: osmo_refdb.annotate_report_levels.tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    def args   = task.ext.args ?: ''
    def levels = split == 'calib' ? '--report-levels-out levels/osmo_refdb.annotate_report_levels.tsv' : ''
    """
    python ${moduleDir}/../templates/benchmark_annotate.py --refs refs --families ${families} --out annotate \\
        --calls-from calls --split ${split} ${levels} ${args}
    for f in calls/*; do
        case "\$f" in *.hmmsearch.tblout|*.blastp.tsv) ;; *) cp -L "\$f" annotate/ ;; esac
    done
    """

    stub:
    """
    mkdir annotate
    if [ "${split}" = calib ]; then mkdir levels && touch levels/osmo_refdb.annotate_report_levels.tsv; fi
    """
}
