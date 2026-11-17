// Simulated read pairs from one split of the references, with per-read truth from CDS overlap.
// split = calib (positive.calib + negative.train: used to decide cascade and report levels)
//       | test  (held-out positives, hard negatives, random background: used to score)
process SIMULATE_READS {
    // calls: wgsim, python + numpy, pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_low'
    tag "${split}"
    input:
    path(refs)
    path(families)
    val(split)

    output:
    path('reads'), emit: reads

    script:
    """
    python ${moduleDir}/templates/09_simulate_reads.py --refs refs --families ${families} \\
        --reads-per-sequence 20 --split ${split} --out reads
    python ${moduleDir}/templates/09b_compute_read_truth.py --reads reads
    """

    stub:
    """
    mkdir reads
    """
}
