//
// The release, the catalog annotated with its HMMs, and the screen that is the gate.
//
include { RELEASE_BUILD          } from '../release_build/main'
include { OSMOTOOL_BUILDCATALOG  } from '../../../modules/local/osmotool/buildcatalog/main'
include { CATALOG_PACKAGE        } from '../../../modules/local/catalog/package/main'
include { CATALOG_SCREEN         } from '../../../modules/local/catalog/screen/main'

workflow RELEASE_CATALOG {
    take:
    families    // families.yaml
    genes       // the catalog's gene file

    main:
    RELEASE_BUILD(families)
    db = RELEASE_BUILD.out.db

    // annotate and index the catalog (osmotool), package it for the release, screen it. build-catalog
    // needs the HMMs with their cutoffs and the DIAMOND database, not the cascade, so it runs alongside
    // the read benchmark
    OSMOTOOL_BUILDCATALOG(genes.map { f -> tuple([id: params.catalog_name], f) }, RELEASE_BUILD.out.db_annotate)
    CATALOG_PACKAGE(params.release, OSMOTOOL_BUILDCATALOG.out.catalog, genes)
    CATALOG_SCREEN(RELEASE_BUILD.out.hmms, CATALOG_PACKAGE.out.catalog, families)

    emit:
    db            = db
    refs          = RELEASE_BUILD.out.refs
    alignments    = RELEASE_BUILD.out.alignments
    hmms          = RELEASE_BUILD.out.hmms
    release_files = RELEASE_BUILD.out.release_files
    cascade       = RELEASE_BUILD.out.cascade
    catalog       = CATALOG_PACKAGE.out.catalog        // [meta, catalogs/<name>/]
    report        = CATALOG_SCREEN.out.report          // the gate report
    done          = RELEASE_BUILD.out.done.mix(CATALOG_SCREEN.out.screen).collect()
}
