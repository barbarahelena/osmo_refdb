// A family's curated Pfam model (families.yaml `pfam_model`, meta.pfam_model) from InterPro, in place
// of an hmmbuild model: its own GA line kept, NAME and ACC set to the family's name.
process PFAM_FETCH {
    tag "${meta.id}"
    label 'process_single'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'refs' (environment.yml next to this file): python + requests
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:ddb8dc66d39a1fc1' :
        'community.wave.seqera.io/library/cd-hit_python_pyyaml_requests:73160cfec41e18d3' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.hmm"), emit: hmm
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests

    script:
    """
    python ${moduleDir}/templates/fetch_pfam_hmm.py ${meta.pfam_model} ${meta.id} ${meta.id}.hmm
    """

    stub:
    """
    printf 'HMMER3/f\\nNAME  ${meta.id}\\nGA    20.0 20.0;\\n//\\n' > ${meta.id}.hmm
    """
}
