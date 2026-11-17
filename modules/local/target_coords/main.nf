// Genomic coordinates of the curated target genes of one genome, from Prodigal's headers.
// The ORF ids were confirmed once against independent annotation (assets/real_genome_validation/target_orfs.tsv).
process TARGET_COORDS {
    // calls: python 3 (stdlib only)
    // image: osmotool (containers/README.md)
    container 'barbarahelena/osmotool:todo'
    label 'process_single'
    tag "${name}"
    input:
    tuple val(name), path(annotation_files)
    path(targets)

    output:
    tuple val(name), path("${name}.gene_regions.tsv"), emit: regions

    script:
    """
    python ${moduleDir}/templates/real_genome_validation/extract_target_coords.py \\
        --prodigal-faa ${name}.prodigal.faa --targets ${targets} --genome ${name} \\
        --out ${name}.gene_regions.tsv
    """

    stub:
    """
    printf 'family\\tcontig\\tstart\\tend\\n' > ${name}.gene_regions.tsv
    """
}
