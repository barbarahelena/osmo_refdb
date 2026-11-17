// HMMs: hmmbuild (or the curated Pfam model), HMM decoys, GA cutoffs, QC scorecard, one pressed database.
process BUILD_HMMS {
    // calls: hmmer (hmmbuild, hmmsearch, hmmpress), mafft, trimal, cd-hit, python + numpy, pyyaml
    // Seqera container 'hmms' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cd-hit_hmmer_trimal_mafft_pruned:cf68077b56599bf6' :
        'community.wave.seqera.io/library/cd-hit_hmmer_trimal_mafft_pruned:608e5325648e4136' }"
    label 'process_medium'
    input:
    path(refs_in, stageAs: 'refs_in')
    path(aln_in, stageAs: 'aln_in/*')   // the families' alignments, the decoys' to reuse, manifest.tsv
    path(families)

    output:
    path('refs'),             emit: refs
    path('alignments'),       emit: alignments
    path('hmms'),             emit: hmms
    path('qc_scorecard.tsv'), emit: qc
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer
    tuple val("${task.process}"), val('mafft'), eval("mafft --version 2>&1 | grep -oE 'v[0-9.]+' | head -1"), topic: versions, emit: versions_mafft
    tuple val("${task.process}"), val('trimal'), eval("trimal --version 2>&1 | grep -oE 'v[0-9][0-9.a-z]*' | head -1 | sed 's/^v//'"), topic: versions, emit: versions_trimal
    tuple val("${task.process}"), val('cd-hit'), eval("cd-hit -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_cd_hit
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pyyaml'), eval("python3 -c 'import yaml; print(yaml.__version__)'"), topic: versions, emit: versions_pyyaml
    tuple val("${task.process}"), val('requests'), eval("python3 -c 'import requests; print(requests.__version__)'"), topic: versions, emit: versions_requests

    script:
    """
    cp -R refs_in refs
    mkdir alignments && cp -L aln_in/* alignments/
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
