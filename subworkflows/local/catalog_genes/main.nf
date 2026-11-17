//
// The gene catalog: download it (skipped if it is already in place, storeDir) and build its
// full competitive index. Needs no release, so it runs alongside the build.
//
include { FETCH_GMGC    } from '../../../modules/local/fetch_gmgc/main'
include { INDEX_CATALOG } from '../../../modules/local/index_catalog/main'

workflow CATALOG_GENES {
    main:
    genes = FETCH_GMGC().map { _f -> params.catalog_genes }
    INDEX_CATALOG(genes)

    emit:
    genes = genes                 // the gene file's path, once it is in place
    index = INDEX_CATALOG.out     // the full index directory
}
