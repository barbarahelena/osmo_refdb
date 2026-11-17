// The reference sequences of the new release: the source release's refs, minus its train/calib/test
// splits, minus study-tagged sequences, plus UniProt sequences for HMM decoys added to families.yaml since.
process RELEASE_REFS {
    // calls: python + requests, pyyaml; cd-hit (strip_ref_tag.py); network (UniProt)
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_medium'
    input:
    path(families)
    path(source, stageAs: 'source')

    output:
    path('refs'), emit: refs

    script:
    """
    mkdir refs
    cp -R source/refs/. refs/
    rm -f refs/*.positive.train.faa refs/*.positive.calib.faa refs/*.positive.test.faa \\
          refs/*.negative.train.faa refs/*.negative.test.faa refs/*.decoy.faa refs/extra_positives_manifest.tsv
    python ${moduleDir}/templates/strip_ref_tag.py --refs refs --tag _study --threads ${task.cpus}
    python ${moduleDir}/templates/01_fetch_refs.py --families ${families} --out refs \\
        --only-hmm-decoys --cache-dir cache
    """

    stub:
    """
    mkdir refs
    """
}
