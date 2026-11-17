//
// The canonical run: everything that makes a release, in dependency order.
//
include { PROVENANCE         } from '../modules/local/provenance/main'
include { CATALOG_GENES      } from '../subworkflows/local/catalog_genes/main'
include { BACDIVE_DATA       } from '../subworkflows/local/bacdive_data/main'
include { RELEASE_CATALOG    } from '../subworkflows/local/release_catalog/main'
include { BACDIVE_ANALYSIS   } from '../subworkflows/local/bacdive_analysis/main'
include { REAL_GENOME_VALIDATION } from '../subworkflows/local/real_genome_validation/main'

workflow OSMO_REFDB {
    main:
    // independent of the release: start at once, run alongside the build
    CATALOG_GENES()
    BACDIVE_DATA()

    // the release, then the catalog on it, then the gate
    RELEASE_CATALOG(CATALOG_GENES.out.genes)
    RELEASE_CATALOG.out.report.view { f -> f.text }

    // after the gate: BacDive on the frozen release, then the NaCl model
    analyzed = channel.empty()
    if( params.gate_passed ) {
        BACDIVE_ANALYSIS(RELEASE_CATALOG.out.release, BACDIVE_DATA.out.proteins)
        REAL_GENOME_VALIDATION(RELEASE_CATALOG.out.release)
        analyzed = BACDIVE_ANALYSIS.out.analyzed.mix(REAL_GENOME_VALIDATION.out.report)
    }

    // the receipt, once everything this invocation does to the release is done
    PROVENANCE(RELEASE_CATALOG.out.release, RELEASE_CATALOG.out.report.mix(analyzed).collect())

    emit:
    release = PROVENANCE.out
}
