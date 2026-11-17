// Read-level recall and precision on one genome's simulated reads: the coordinate-based truth
// (from wgsim's read names and the target genes' regions) against DIAMOND (osmotool profile's
// alignments) and HMM (hmmscan on the reads' ORFs), per target family.
process VALIDATION_READRECALL {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(reads), path(regions), path(blastx), path(tblout)

    output:
    tuple val(meta), path("${meta.id}.read_truth.tsv"), path("${meta.id}.real_recall_comparison.tsv"), emit: results
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    def T = "${moduleDir}/templates/real_genome_validation"
    """
    python ${T}/real_read_truth.py \\
        --r1 ${reads[0]} --r2 ${reads[1]} --regions ${regions} --read-length 150 --out ${meta.id}.read_truth.tsv
    python ${T}/compare_real_recall.py \\
        --truth ${meta.id}.read_truth.tsv --diamond-blastx ${blastx} \\
        --hmm-tblout ${tblout} --out ${meta.id}.real_recall_comparison.tsv
    """

    stub:
    """
    touch ${meta.id}.read_truth.tsv ${meta.id}.real_recall_comparison.tsv
    """
}
