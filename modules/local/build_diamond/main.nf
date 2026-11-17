// DIAMOND database from the train references plus decoy and fused-ORF references, its per-family
// bit-score cutoffs, and the lists of families kept out of reported output.
process BUILD_DIAMOND {
    // calls: diamond, python + numpy, pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_medium'
    input:
    path(refs_in, stageAs: 'refs_in')
    path(families)

    output:
    path('refs'),          emit: refs
    path('release_files'), emit: release_files

    script:
    """
    cp -R refs_in refs
    mkdir release_files
    T=${moduleDir}/templates
    python \$T/08a_build_decoy_refs.py --refs refs --families ${families}
    python \$T/08d_build_fusion_refs.py --refs refs --families ${families}
    THREADS=${task.cpus} bash \$T/08_build_diamond_db.sh refs release_files osmo_refdb
    python \$T/08b_calibrate_diamond_cutoffs.py --refs refs --release release_files \\
        --release-name osmo_refdb --families ${families} --threads ${task.cpus}
    python \$T/08c_write_scope_manifest.py --families ${families} \\
        --profile-out release_files/osmo_refdb.profile_excluded_families.txt \\
        --annotate-out release_files/osmo_refdb.annotate_excluded_families.txt
    """

    stub:
    """
    cp -R refs_in refs
    mkdir release_files && touch release_files/osmo_refdb.dmnd
    """
}
