// The families kept out of osmotool's reported output, as the flat lists profile and annotate read
// (--exclude_families): families with scope: annotate_only (out of profile only: murB is in the panel
// for annotate's co-occurrence checks, meaningless in read profiling), and every decoy, DIAMOND
// (<family>_decoy, decoy_from_negatives) and HMM (hmm_decoys), out of both. They stay in the
// databases: a decoy exists to win a protein away from a wrong call, not to be counted.
process DIAMOND_SCOPE {
    tag "${name}"
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    val(name)
    path(families)

    output:
    path("${name}.{profile,annotate}_excluded_families.txt"), emit: excluded
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    python3 <<'PYEOF'
    import yaml

    with open("${families}") as fh:
        data = yaml.safe_load(fh)
    annotate_only = [fam["name"] for fam in data["families"] if fam.get("scope") == "annotate_only"]
    decoys = [f"{fam['name']}_decoy" for fam in data["families"] if fam.get("decoy_from_negatives")]
    decoys += [d["name"] for d in data.get("hmm_decoys") or []]

    for mode, excluded in (("profile", annotate_only + decoys), ("annotate", decoys)):
        with open(f"${name}.{mode}_excluded_families.txt", "w") as fh:
            fh.write("".join(f"{x}\\n" for x in excluded))
        print(f"{mode}-excluded: {len(excluded)} families {excluded}")
    PYEOF
    """

    stub:
    """
    touch ${name}.profile_excluded_families.txt ${name}.annotate_excluded_families.txt
    """
}
