//
// BacDive strains and genomes, then Prodigal proteins and the random-panel pool search for every
// genome (an array of nchunks tasks). Release-independent: runs alongside the build.
//
include { FETCH_BACDIVE    } from '../../../modules/local/fetch_bacdive/main'
include { BACDIVE_PROTEINS } from '../../../modules/local/bacdive_proteins/main'

workflow BACDIVE_DATA {
    main:
    data     = FETCH_BACDIVE()
    chunks   = channel.of(0..(params.nchunks - 1))
    proteins = BACDIVE_PROTEINS(chunks, data).collect()

    emit:
    proteins = proteins           // one value, once every chunk is done
}
