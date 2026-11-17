// The release's alignments/ and the HMM decoys' input sets, as the release keeps them: the alignments
// made in this build, those reused from the source release (which are not any step's output, so they
// pass through here to be published), and manifest.tsv, what each one is (the next build reads it to
// decide what it can reuse). Published to releases/<release>/{alignments,refs} (conf/modules.config).
process ALIGNMENTS_RECORD {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    path(alignments, stageAs: 'in_alignments/*')   // families' and decoys' alignments, manifest.tsv
    path(decoy_refs, stageAs: 'in_refs/*')         // the decoys' input and clustered sets

    output:
    path('alignments'), emit: alignments
    path('refs/*'),     optional: true, emit: decoy_refs

    script:
    """
    mkdir alignments refs
    cp -L in_alignments/* alignments/
    if [ -d in_refs ]; then cp -L in_refs/* refs/; fi
    """

    stub:
    """
    mkdir alignments refs
    cp -L in_alignments/* alignments/
    if [ -d in_refs ]; then cp -L in_refs/* refs/; fi
    """
}
