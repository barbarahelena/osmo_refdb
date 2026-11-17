// Make releases/<release>/phenotype_validation point at the release-independent data
// (genomes, cached proteins, random panel), after checking that every genome has proteins.
process PREPARE_BACDIVE {
    // calls: bash, ls, grep, ln, cp
    // image: any small image (python + pyyaml as for SCREEN will do)
    container 'todo.invalid/prepare_bacdive'
    label 'process_single'
    input:
    val(release)
    val(proteins_done)

    output:
    val(release)

    script:
    """
    cd ${params.repo}
    N_GEN=\$(ls ${params.data_dir}/genomes | grep -c '\\.fna' || true)
    N_PROT=\$(ls ${params.data_dir}/random_panel/strains/*/*.pool_counts.tsv 2>/dev/null | wc -l)
    if [ "\$N_PROT" -lt "\$N_GEN" ]; then
        echo "ERROR: proteins for \$N_PROT of \$N_GEN genomes" >&2
        exit 1
    fi
    bash hpc/prepare_bacdive_release.sh ${release} ${params.data_dir}
    """

    stub:
    """
    echo stub
    """
}
