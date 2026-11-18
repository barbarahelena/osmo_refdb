// The sequences the DIAMOND database is built from, <name>.train_refs.faa (it ships with the release:
// osmotool's marker normalisation reads it): every family's train positives without the
// _refseq-tagged ones, then the decoy references, then the fused-ORF references, each in name order.
process DIAMOND_TRAINREFS {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.9--he361c42_0' :
        'biocontainers/diamond:2.2.9--he361c42_0' }"

    input:
    val(name)
    path(train, stageAs: 'train/*')     // SPLIT_REFS: <family>.positive.train.faa
    path(extra, stageAs: 'extra/*')     // DIAMOND_EXTRAREFS: *.decoy.faa, *_fused.faa

    output:
    path("${name}.train_refs.faa"), emit: fasta

    script:
    """
    export LC_ALL=C   # name order, whatever the locale
    OUT=${name}.train_refs.faa
    : > \$OUT
    for f in train/*.positive.train.faa; do
        [ -s "\$f" ] || continue
        awk '/^>/ { keep = (\$0 !~ /_refseq\\|/) } keep' "\$f" >> \$OUT
    done
    for f in extra/*.decoy.faa extra/*_fused.faa; do
        if [ -e "\$f" ]; then cat "\$f" >> \$OUT; fi
    done
    """

    stub:
    """
    printf '>a\\nMK\\n' > ${name}.train_refs.faa
    """
}
