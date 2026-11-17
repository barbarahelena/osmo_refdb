//
// After the gate: BacDive annotation with the frozen release's HMMs (an array of nchunks tasks),
// the analysis, and the NaCl model.
//
include { PREPARE_BACDIVE  } from '../../../modules/local/prepare_bacdive/main'
include { BACDIVE_ANNOTATE } from '../../../modules/local/bacdive_annotate/main'
include { ANALYZE_NACL     } from '../../../modules/local/analyze_nacl/main'

workflow BACDIVE_ANALYSIS {
    take:
    release     // the release's name, once it is final
    proteins    // once every genome has its proteins

    main:
    prepared  = PREPARE_BACDIVE(release, proteins)
    chunks    = channel.of(0..(params.nchunks - 1))
    annotated = BACDIVE_ANNOTATE(chunks, prepared).collect()
    analyzed  = ANALYZE_NACL(prepared, annotated)

    emit:
    analyzed = analyzed
}
