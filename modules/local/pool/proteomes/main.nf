// The screening proteomes for the random pool: Prodigal's proteins of each genome (-p meta, as the
// BacDive steps call them), headed <bacdive_id>|<protein id>, in one FASTA.
process POOL_PROTEOMES {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prodigal:2.6.3--h577a1d6_11' :
        'biocontainers/prodigal:2.6.3--h577a1d6_11' }"

    input:
    path(genomes)       // BACDIVE_GENOMES: <bid>_<accession>.fna[.gz]

    output:
    path('screen_proteomes.faa'), emit: proteomes
    tuple val("${task.process}"), val('prodigal'), eval("prodigal -v 2>&1 | grep -oE 'V[0-9.]+' | tr -d V"), topic: versions, emit: versions_prodigal

    script:
    """
    mkdir proteins
    call() {   # <genome>: proteins/<bid>.faa
        b=\$(basename "\$1"); bid=\${b%%_*}
        zcat -f "\$1" > "proteins/\$bid.fna"
        prodigal -i "proteins/\$bid.fna" -a "proteins/\$bid.raw.faa" -o /dev/null -p meta -q
        sed "s/^>/>\$bid|/" "proteins/\$bid.raw.faa" > "proteins/\$bid.faa"
        rm "proteins/\$bid.fna" "proteins/\$bid.raw.faa"
    }
    export -f call
    ls ${genomes}/*.fna* | xargs -P ${task.cpus} -I{} bash -c 'call "\$1"' _ {}
    cat \$(ls proteins/*.faa | sort) > screen_proteomes.faa
    """

    stub:
    """
    printf '>1|c_1\\nMK\\n' > screen_proteomes.faa
    """
}
