// osmotool annotate of one reference genome against the release (DIAMOND + HMM, Prodigal proteins kept).
// Used for the system calls on the reference genomes and as the first step of the read-recall check.
process ANNOTATE_GENOME {
    // calls: osmotool annotate (prodigal, hmmsearch, diamond)
    // image: osmotool (containers/README.md)
    container 'barbarahelena/osmotool:todo'
    label 'process_low'
    tag "${name}"
    input:
    val(release)
    tuple val(name), path(genome)

    output:
    tuple val(name), path("${name}.*"), emit: annotated

    script:
    """
    osmotool annotate ${params.repo}/releases/${release} ${genome} \\
        --out_prefix ${name} --threads ${task.cpus} --keep_proteins
    """

    stub:
    """
    printf 'system\\tpresent\\tsubunit_count\\trescued_hits\\n' > ${name}.systems.tsv
    printf 'gene\\traw_count\\n' > ${name}.gene_counts.tsv
    echo ">${name}_1 # 1 # 300 # 1 # stub" > ${name}.prodigal.faa
    """
}
