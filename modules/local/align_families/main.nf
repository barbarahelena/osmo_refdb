// The alignments the HMMs are built from. The source release's alignments are reused (multi-threaded
// MAFFT is not reproducible run to run); only families whose train set differs from the source's are realigned.
process ALIGN_FAMILIES {
    // calls: mafft, trimal, python + pyyaml
    // Seqera container 'align' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/trimal_mafft_python_pyyaml:0aafcb62c9cac1a6' :
        'community.wave.seqera.io/library/trimal_mafft_python_pyyaml:d3ab7e2b4cc1b639' }"
    label 'process_medium'
    input:
    path(refs)
    path(families)
    path(source, stageAs: 'source')

    output:
    path('alignments'), emit: alignments
    tuple val("${task.process}"), val('mafft'), eval("mafft --version 2>&1 | grep -oE 'v[0-9.]+' | head -1"), topic: versions, emit: versions_mafft
    tuple val("${task.process}"), val('trimal'), eval("trimal --version 2>&1 | grep -oE 'v[0-9][0-9.a-z]*' | head -1 | sed 's/^v//'"), topic: versions, emit: versions_trimal
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    mkdir alignments
    cp source/alignments/* alignments/
    CHANGED=""
    for f in refs/*.positive.train.faa; do
        cmp -s \$f source/refs/\$(basename \$f) || CHANGED="\$CHANGED \$(basename \$f .positive.train.faa)"
    done
    if [ -n "\$CHANGED" ]; then
        echo "Train set changed, realigning:\$CHANGED"
        ALIGN_ONLY="\${CHANGED# }" THREADS=${task.cpus} bash ${moduleDir}/templates/04_align_trim.sh refs alignments ${families}
    fi
    """

    stub:
    """
    mkdir alignments
    """
}
