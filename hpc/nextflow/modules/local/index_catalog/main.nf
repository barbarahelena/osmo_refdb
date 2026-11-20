// The full, competitive index of the catalog (classic bwa: bwa-mem3 would need 100+ GB).
// It depends on the gene file and the aligner only, so it runs alongside the build.
// storeDir: built here, moved to full_index_dir when complete (never half there); skipped if it exists.
process INDEX_CATALOG {
    // calls: osmotool index-catalog --aligner bwa (classic bwa)
    // image: the osmotool image, which must also contain classic bwa
    container 'todo.invalid/index_catalog'
    label 'process_single'
    storeDir "${file(params.full_index_dir).parent}"
    input:
    val(genes)

    output:
    path("${file(params.full_index_dir).name}")

    script:
    """
    osmotool index-catalog --aligner bwa --genes ${genes} -o ${file(params.full_index_dir).name}
    test -f ${file(params.full_index_dir).name}/catalog_index.json
    """

    stub:
    """
    mkdir ${file(params.full_index_dir).name}; echo stub > ${file(params.full_index_dir).name}/catalog_index.json
    """
}
