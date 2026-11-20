// Pfam-A, pinned to release 39.0 (the release the random-panel pool was drawn from), for the pool's
// candidate families. Pfam.version records it.
process POOL_PFAM {
    label 'process_single'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // Seqera container 'bacdive-fetch' (environment.yml next to this file): curl
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/curl_findutils_gawk_gzip_pruned:dfd51b437477d4e7' :
        'community.wave.seqera.io/library/curl_findutils_gawk_gzip_pruned:47915284675a72c3' }"

    input:
    val(release)    // the Pfam release, e.g. 39.0

    output:
    path('Pfam-A.hmm.gz'), emit: pfam
    path('Pfam.version'),  emit: version

    script:
    """
    BASE=https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam${release}
    curl -sS -f -L --retry 3 -o Pfam-A.hmm.gz \$BASE/Pfam-A.hmm.gz
    curl -sS -f -L --retry 3 -o Pfam.version.gz \$BASE/Pfam.version.gz
    python3 -c "import gzip; open('Pfam.version', 'wb').write(gzip.open('Pfam.version.gz').read())"
    cat Pfam.version
    """

    stub:
    """
    touch Pfam-A.hmm.gz
    echo "Pfam release ${release}" > Pfam.version
    """
}
