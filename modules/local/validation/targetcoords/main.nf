// Genomic coordinates of the curated target genes of one genome, from Prodigal's headers.
// The ORF ids were confirmed once against independent annotation (assets/real_genome_validation/target_orfs.tsv).
process VALIDATION_TARGETCOORDS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(proteins)   // osmotool annotate --keep_proteins: <genome>.prodigal.faa
    path(targets)

    output:
    tuple val(meta), path("${meta.id}.gene_regions.tsv"), emit: regions
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python ${moduleDir}/templates/real_genome_validation/extract_target_coords.py \\
        --prodigal-faa ${proteins} --targets ${targets} --genome ${meta.id} \\
        --out ${meta.id}.gene_regions.tsv
    """

    stub:
    """
    printf 'family\\tcontig\\tstart\\tend\\n' > ${meta.id}.gene_regions.tsv
    """
}
