// The DIAMOND database's references besides the train sets: decoy references from the negatives of
// families marked decoy_from_negatives (templates/08a_build_decoy_refs.py), and fused-ORF references
// for the fusion_partner pairs (templates/08d_build_fusion_refs.py). Published into the release's refs/.
process DIAMOND_EXTRAREFS {
    label 'process_single'

    // Seqera container 'py-yaml' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pyyaml:9b2a7bd921471580' :
        'community.wave.seqera.io/library/python_pyyaml:9e6e168c801c7a81' }"

    input:
    path(refs)          // SPLIT_REFS
    path(families)

    output:
    path('extra/*.faa'), optional: true, emit: refs   // <family>.decoy.faa, <a>_<b>_fused.faa
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    python ${moduleDir}/templates/08a_build_decoy_refs.py --refs ${refs} --families ${families} --out extra
    python ${moduleDir}/templates/08d_build_fusion_refs.py --refs ${refs} --families ${families} --out extra
    """

    stub:
    """
    mkdir extra
    printf '>betL|decoy\\nMK\\n' > extra/betL.decoy.faa
    """
}
