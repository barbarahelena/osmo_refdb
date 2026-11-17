// osmotool annotate with the release's HMMs, per genome, from the cached proteins. Links the release's
// phenotype_validation to the shared data first (see templates/bacdive_annotate.sh).
process BACDIVE_ANNOTATE {
    // calls: osmotool annotate (hmmsearch, diamond, prodigal)
    // image: osmotool (containers/README.md)
    container 'barbarahelena/osmotool:todo'
    label 'process_low'
    tag "chunk ${chunk}/${params.nchunks}"
    input:
    val(chunk)
    val(release)
    val(proteins_done)

    output:
    val(chunk)

    script:
    """
    cd ${params.repo}
    # Point the release's phenotype_validation at the release-independent data
    # (genomes, cached proteins, random panel). Idempotent: every chunk does it.
    DATA=\$(cd ${params.data_dir} && pwd)
    PV=releases/${release}/phenotype_validation
    mkdir -p \$PV/strategy_control
    for d in genomes random_panel; do
        ln -sfn \$DATA/\$d \$PV/\$d
        if [ -e \$DATA/strategy_control/\$d ]; then ln -sfn \$DATA/strategy_control/\$d \$PV/strategy_control/\$d; fi
    done
    for f in pilot_strains.tsv bacdive_classes.tsv; do
        if [ -s \$DATA/\$f ]; then cp -f \$DATA/\$f \$PV/\$f; fi
    done
    if [ -s \$DATA/strategy_control/pilot_strains.tsv ]; then
        cp -f \$DATA/strategy_control/pilot_strains.tsv \$PV/strategy_control/
    fi

    # osmotool annotate per genome, from the cached proteins; task i takes every nchunks-th genome
    for D in \$PV \$PV/strategy_control; do
        [ -d \$D/genomes ] || continue
        CHUNK=${chunk}/${params.nchunks} bash ${moduleDir}/templates/phenotype_validation/run_annotate_cached.sh \$D releases/${release} ${task.cpus}
    done
    """

    stub:
    """
    echo stub
    """
}
