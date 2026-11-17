//
// BacDive strains and genomes, then Prodigal proteins and the random-panel pool search for every
// genome (nchunks tasks, each every nchunks-th genome). Release-independent: runs alongside the build.
//
include { BACDIVE_FETCH    } from '../../../modules/local/bacdive/fetch/main'
include { BACDIVE_PROTEINS } from '../../../modules/local/bacdive/proteins/main'

workflow BACDIVE_DATA {
    take:
    snapshot        // harvested BacDive batches to start from (bacdive_raw/), or []
    pool            // the random-panel pool, pool.hmm

    main:
    BACDIVE_FETCH(snapshot)
    chunks = channel.of(0..(params.nchunks - 1))
    BACDIVE_PROTEINS(chunks, BACDIVE_FETCH.out.genomes, BACDIVE_FETCH.out.control_genomes, pool)

    emit:
    genomes         = BACDIVE_FETCH.out.genomes
    control_genomes = BACDIVE_FETCH.out.control_genomes
    pilot           = BACDIVE_FETCH.out.pilot
    classes         = BACDIVE_FETCH.out.classes
    control_pilot   = BACDIVE_FETCH.out.control_pilot
    chunks          = BACDIVE_PROTEINS.out.chunk                         // [chunk, strains/, control_strains/]
    strains         = BACDIVE_PROTEINS.out.strains.flatten().toList()    // every strain's <bid>/
    control_strains = BACDIVE_PROTEINS.out.control_strains.flatten().toList()
}
