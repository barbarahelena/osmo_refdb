// Paired reads straight from a real genome with wgsim. The read headers carry the true coordinates,
// which is what makes the ground truth independent of the database (real_read_truth.py).
process SIMULATE_GENOME_READS {
    // calls: wgsim
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_single'
    tag "${name}"
    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path("${name}_R1.fastq"), path("${name}_R2.fastq"), emit: reads

    script:
    """
    wgsim -S 1 -N ${params.recall_read_pairs} -1 150 -2 150 ${genome} ${name}_R1.fastq ${name}_R2.fastq > /dev/null
    """

    stub:
    """
    touch ${name}_R1.fastq ${name}_R2.fastq
    """
}
