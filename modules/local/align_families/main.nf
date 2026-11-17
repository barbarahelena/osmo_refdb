// The alignments the HMMs are built from. The source release's alignments are reused (multi-threaded
// MAFFT is not reproducible run to run); only families whose train set differs from the source's are realigned.
process ALIGN_FAMILIES {
    // calls: mafft, trimal, python + pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_medium'
    input:
    path(refs)
    path(families)
    path(source, stageAs: 'source')

    output:
    path('alignments'), emit: alignments

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
