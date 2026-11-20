// The current UniProt release (the X-UniProt-Release header of a size=0 search). It keys the stored
// fetch (--refs_store): within one release a query's population is fixed, so the same families file
// fetches the same references. Never cached: UniProt can release between two runs.
process UNIPROT_RELEASE {
    label 'process_single'
    cache false
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'refs' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:bb86db01ccb6fea3' :
        'community.wave.seqera.io/library/cd-hit_diamond_procps-ng_python_pruned:8b81aff3d2d73d3c' }"

    output:
    stdout emit: release
    // the 'refs' image's versions for REFS_FETCHFAMILY and REFS_FETCHSHARED too, which cannot emit them (storeDir)
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests

    script:
    """
    python3 - <<'PYEOF'
    import sys, requests
    r = requests.get("https://rest.uniprot.org/uniprotkb/search", params={"query": "reviewed:true", "size": 0}, timeout=30)
    r.raise_for_status()
    release = r.headers.get("X-UniProt-Release")
    if not release:
        sys.exit("ERROR: no X-UniProt-Release header")
    print(release, end="")
    PYEOF
    """

    stub:
    """
    printf 'stub_release'
    """
}
