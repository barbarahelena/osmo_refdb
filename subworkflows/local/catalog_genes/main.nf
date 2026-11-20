//
// The gene catalog: download it (skipped if it is already in place, storeDir). Needs no release, so it
// runs alongside the build; RELEASE_CATALOG annotates and indexes it for each release.
//
include { CATALOG_FETCH } from '../../../modules/local/catalog/fetch/main'

workflow CATALOG_GENES {
    main:
    CATALOG_FETCH()

    emit:
    genes = CATALOG_FETCH.out.genes   // the gene file, once it is in place
}
