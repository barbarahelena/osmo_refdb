// The release's hmms/: each family's GA cutoff from its calibration scores, and each decoy's fixed GA
// (templates/calibrate_cutoffs.py, which writes the GA lines and cutoff_manifest.tsv), the QC scorecard
// (templates/qc_scorecard.py), and every model in one pressed database, hmms/osmo_refdb.hmm.
process HMM_CALIBRATE {
    label 'process_single'

    // Seqera container 'hmms' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_hmmer_trimal_mafft_pruned:cf68077b56599bf6' :
        'community.wave.seqera.io/library/cd-hit_hmmer_trimal_mafft_pruned:608e5325648e4136' }"

    input:
    path(models, stageAs: 'models/*')     // every family's and decoy's .hmm (HMMER_HMMBUILD, PFAM_FETCH)
    path(scores, stageAs: 'scores/*')     // HMMER_HMMSEARCH: <family>[.<decoy>].{positive,negative}.tblout
    path(refs)                            // SPLIT_REFS: the purity and length-outlier manifests
    path(families)

    output:
    path('hmms'),             emit: hmms
    path('qc_scorecard.tsv'), emit: qc
    tuple val("${task.process}"), val('hmmer'), eval("hmmpress -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml

    script:
    """
    # the GA lines are written into the models: copies, not the other steps' files
    mkdir -p hmms/scores
    cp -L models/*.hmm hmms/
    cp -L scores/* hmms/scores/
    python ${moduleDir}/templates/calibrate_cutoffs.py --hmms hmms --families ${families}
    python ${moduleDir}/templates/qc_scorecard.py --refs ${refs} --hmms hmms --families ${families} --out qc_scorecard.tsv
    # one database of every model, in name order
    cat hmms/*.hmm > osmo_refdb.hmm
    mv osmo_refdb.hmm hmms/osmo_refdb.hmm
    hmmpress hmms/osmo_refdb.hmm > /dev/null
    """

    stub:
    """
    mkdir -p hmms/scores
    cp -L models/*.hmm hmms/
    touch hmms/osmo_refdb.hmm hmms/osmo_refdb.hmm.h3m hmms/cutoff_manifest.tsv qc_scorecard.tsv
    """
}
