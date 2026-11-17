//
// The release, the catalog annotated with its HMMs, and the screen that is the gate.
//
include { RELEASE_BUILD    } from '../release_build/main'
include { ANNOTATE_CATALOG } from '../../../modules/local/annotate_catalog/main'
include { SCREEN           } from '../../../modules/local/screen/main'

workflow RELEASE_CATALOG {
    take:
    genes     // the gene file's path, once it is in place

    main:
    RELEASE_BUILD()
    cat     = ANNOTATE_CATALOG(RELEASE_BUILD.out.release, genes)
    SCREEN(cat)

    emit:
    release = cat                 // the release's name, once the catalog is in it
    report  = SCREEN.out.report   // the gate report
}
