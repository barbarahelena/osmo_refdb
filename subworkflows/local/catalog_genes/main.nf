//
// The gene catalog: download it (skipped if it is already in place, storeDir) and build its
// full competitive index (storeDir too, conf/modules.config). Needs no release, so it runs
// alongside the build.
//
include { CATALOG_FETCH         } from '../../../modules/local/catalog/fetch/main'
include { OSMOTOOL_INDEXCATALOG } from '../../../modules/local/osmotool/indexcatalog/main'

workflow CATALOG_GENES {
    main:
    CATALOG_FETCH()
    genes = CATALOG_FETCH.out.genes
    OSMOTOOL_INDEXCATALOG(genes.map { f -> tuple([id: params.catalog_name], f) })

    emit:
    genes    = genes                              // the gene file, once it is in place
    index    = OSMOTOOL_INDEXCATALOG.out.index    // the full index directory
}
