// What this run is built from: the code (commit and a checksum of any uncommitted
// change). Looked at again on every run (cache false) and handed to the release build, so
// changing the code changes that task's inputs and -resume does not reuse a build made from
// something else. (A changed image changes the task's hash by itself: the container is part of it.)
// Runs on the node itself: git and sha256sum are on every cluster.
process PIN {
    label 'process_single'
    cache false
    output:
    tuple env('COMMIT'), env('DIRTY'), env('DIFF')

    script:
    """
    COMMIT=\$(git -C ${params.repo} rev-parse HEAD)
    DIRTY=\$(git -C ${params.repo} status --porcelain --untracked-files=no | wc -l | tr -d ' ')
    DIFF=\$(git -C ${params.repo} diff HEAD | sha256sum | cut -c1-16)
    """

    stub:
    """
    COMMIT=stubcommit; DIRTY=0; DIFF=stubdiff
    """
}
