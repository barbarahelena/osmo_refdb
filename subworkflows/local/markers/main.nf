//
// The release's markers/: universal single-copy marker genes (markers.yaml), for osmotool profile's
// copies-per-genome normalisation. Built apart from the osmo families, so they never compete with them
// in a best-hit contest: their own references, split and DIAMOND database, then packaged.
//   REFS_BUILD (no random background: no benchmark) -> SPLIT_REFS -> DIAMOND_BUILD -> MARKERS_PACKAGE
// The marker manifest's DIAMOND precision/recall columns stay empty: the markers have no benchmark here.
//
include { REFS_BUILD      } from '../refs_build/main'
include { SPLIT_REFS      } from '../../../modules/local/split_refs/main'
include { DIAMOND_BUILD   } from '../diamond_build/main'
include { MARKERS_PACKAGE } from '../../../modules/local/markers/package/main'

workflow MARKERS {
    take:
    markers     // markers.yaml

    main:
    REFS_BUILD(markers, false)
    SPLIT_REFS(REFS_BUILD.out.refs, markers)
    DIAMOND_BUILD(markers, SPLIT_REFS.out.refs, SPLIT_REFS.out.train, SPLIT_REFS.out.faa)
    MARKERS_PACKAGE(DIAMOND_BUILD.out.release_files, REFS_BUILD.out.manifest, markers)

    emit:
    markers = MARKERS_PACKAGE.out.markers   // the release's markers/
}
