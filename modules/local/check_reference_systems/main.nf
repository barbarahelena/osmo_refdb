// Compare the system calls on the reference genomes with the committed expectations
// (assets/real_genome_validation/expected_systems.tsv). A failed expectation is reported, not fatal:
// the verdict is part of the release's validation, and the workflow goes on.
process CHECK_REFERENCE_SYSTEMS {
    // calls: python 3 (stdlib only)
    // image: osmotool (containers/README.md)
    container 'barbarahelena/osmotool:todo'
    label 'process_single'
    input:
    path(annotation_files)
    path(expected)

    output:
    path('reference_systems_check.*'), emit: report

    script:
    """
    # the checker reads <release>/real_genome_validation/reference_genomes/<genome>/<genome>.*.tsv
    D=release/real_genome_validation/reference_genomes
    for f in *.systems.tsv *.gene_counts.tsv; do
        g=\${f%%.*}
        mkdir -p \$D/\$g
        cp \$f \$D/\$g/
    done
    set +e
    python ${moduleDir}/templates/real_genome_validation/check_reference_systems.py --release release \\
        --expected ${expected} \\
        > reference_systems_check.txt
    echo "exit status \$? (1 = an expectation failed)" >> reference_systems_check.txt
    cp release/real_genome_validation/reference_systems_check.tsv .
    """

    stub:
    """
    echo stub > reference_systems_check.txt
    """
}
