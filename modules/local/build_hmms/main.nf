// HMMs: hmmbuild (or the curated Pfam model), HMM decoys, GA cutoffs, QC scorecard, one pressed database.
process BUILD_HMMS {
    // calls: hmmer (hmmbuild, hmmsearch, hmmpress), mafft, trimal, cd-hit, python + numpy, pyyaml
    // image: osmo-build (containers/README.md)
    container 'barbarahelena/osmo-build:todo'
    label 'process_medium'
    input:
    path(refs_in, stageAs: 'refs_in')
    path(aln_in, stageAs: 'aln_in')
    path(families)

    output:
    path('refs'),             emit: refs
    path('alignments'),       emit: alignments
    path('hmms'),             emit: hmms
    path('qc_scorecard.tsv'), emit: qc

    script:
    """
    cp -R refs_in refs
    mkdir alignments && cp -R aln_in/. alignments/
    T=${moduleDir}/templates
    THREADS=${task.cpus} bash \$T/05_build_hmms.sh refs alignments hmms ${families}
    python \$T/05b_build_hmm_decoys.py --refs refs --aln alignments --hmms hmms --families ${families} \\
        --threads ${task.cpus} --reuse-alignments
    python \$T/06_calibrate_cutoffs.py --hmms hmms --families ${families}
    python \$T/06b_qc_scorecard.py --refs refs --hmms hmms --families ${families} --out qc_scorecard.tsv
    bash \$T/07_press_hmms.sh hmms osmo_refdb
    """

    stub:
    """
    cp -R refs_in refs
    mkdir alignments hmms && touch hmms/osmo_refdb.hmm qc_scorecard.tsv
    """
}
