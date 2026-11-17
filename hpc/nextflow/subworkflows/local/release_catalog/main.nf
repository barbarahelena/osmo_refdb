//
// The release, the catalog annotated with its HMMs, and the screen that is the gate.
//
include { BUILD_RELEASE    } from '../../../modules/local/build_release/main'
include { ANNOTATE_CATALOG } from '../../../modules/local/annotate_catalog/main'
include { SCREEN           } from '../../../modules/local/catalog_screen/main'

workflow RELEASE_CATALOG {
    take:
    pin       // what the build is made from (PIN)
    genes     // the gene file's path, once it is in place

    main:
    release = BUILD_RELEASE(pin)
    cat     = ANNOTATE_CATALOG(release, genes)
    SCREEN(cat)

    emit:
    release = cat                 // the release's name, once the catalog is in it
    report  = SCREEN.out.report   // the gate report
}
