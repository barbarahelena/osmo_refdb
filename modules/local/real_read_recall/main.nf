// DIAMOND (osmotool profile) and HMM (orfm + hmmscan) on the genome's reads, scored against the
// coordinate-based truth: recall and precision per target family, DIAMOND next to HMM.
process REAL_READ_RECALL {
    // calls: osmotool profile (diamond, hmmsearch, orfm), hmmscan, orfm, python 3
    // image: osmotool (containers/README.md)
    container 'barbarahelena/osmotool:todo'
    label 'process_medium'
    tag "${name}"
    input:
    val(release)
    tuple val(name), path(r1), path(r2), path(regions)

    output:
    tuple val(name), path("${name}.*"), emit: results

    script:
    """
    REL=${params.repo}/releases/${release}
    EXTRA=""
    if osmotool profile --help 2>/dev/null | grep -q -- "--no_markers"; then EXTRA="--no_markers"; fi
    osmotool profile \$REL -1 ${r1} -2 ${r2} --out_prefix ${name} --threads ${task.cpus} \\
        --min_identity 0.40 --min_query_cover 0.50 --keep_aln \$EXTRA
    bash ${moduleDir}/templates/translate_and_hmmscan.sh ${r1} ${r2} \$REL/hmms/osmo_refdb.hmm ${name}.reads ${task.cpus}
    python ${moduleDir}/templates/real_genome_validation/real_read_truth.py \\
        --r1 ${r1} --r2 ${r2} --regions ${regions} --read-length 150 --out ${name}.read_truth.tsv
    python ${moduleDir}/templates/real_genome_validation/compare_real_recall.py \\
        --truth ${name}.read_truth.tsv --diamond-blastx ${name}.blastx.tsv \\
        --hmm-tblout ${name}.reads.hmmscan.tblout --out ${name}.real_recall_comparison.tsv
    """

    stub:
    """
    echo stub > ${name}.real_recall_comparison.tsv
    """
}
