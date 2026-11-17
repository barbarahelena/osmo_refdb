// osmotool annotate on held-out full-length proteins (HMM and DIAMOND), scored per family.
process BENCHMARK_ANNOTATE {
    // calls: osmotool annotate (prodigal, hmmsearch, diamond), python + pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_medium'
    input:
    path(refs)
    path(hmms)
    path(release_files)
    path(cascade)
    path(families)

    output:
    path('annotate'), emit: annotate

    script:
    """
    mkdir rel
    for f in release_files/* cascade/*; do [ -e \$f ] && ln -s \$PWD/\$f rel/\$(basename \$f); done
    ln -s \$PWD/hmms rel/hmms
    python ${moduleDir}/templates/11b_benchmark_annotate.py --refs refs --release rel \\
        --families ${families} --out annotate --threads ${task.cpus}
    """

    stub:
    """
    mkdir annotate
    """
}
