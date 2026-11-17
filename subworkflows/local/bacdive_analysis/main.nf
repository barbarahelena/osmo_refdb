//
// After the gate: BacDive annotation with the frozen release's HMMs (an array of nchunks tasks),
// the analysis, and the NaCl model.
//
include { BACDIVE_ANNOTATE } from '../../../modules/local/bacdive_annotate/main'
include { ANALYZE_NACL     } from '../../../modules/local/analyze_nacl/main'

workflow BACDIVE_ANALYSIS {
    take:
    release     // the release's name, once it is final
    proteins    // once every genome has its proteins (checked in BACDIVE_DATA)

    main:
    chunks    = channel.of(0..(params.nchunks - 1))
    annotated = BACDIVE_ANNOTATE(chunks, release, proteins).collect()
    analyzed  = ANALYZE_NACL(release, annotated)

    emit:
    analyzed = analyzed
}
