// The provenance receipt, written last, so it covers the catalog (and the BacDive and NaCl results,
// when the gate has been passed) as well as the build. It probes the version of every tool, so it runs in
// the image that has all of them.
process PROVENANCE {
    // calls: python, git, and every tool's --version (write_provenance.py probes them)
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_single'
    input:
    val(release)
    val(after)

    output:
    val(release)

    script:
    """
    cd ${params.repo}
    python ${moduleDir}/templates/write_provenance.py --release releases/${release} \\
        --repo ${params.repo} --families ${params.families_file}
    test -s releases/${release}/PROVENANCE.json
    """

    stub:
    """
    echo stub
    """
}
