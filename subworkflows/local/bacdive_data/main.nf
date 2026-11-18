//
// The BacDive strains, their genomes, and per genome its Prodigal proteins and random-pool search.
// Release-independent: runs alongside the build.
//   default         the 43 pilot strains (assets/bacdive/pilot_strains.tsv): a quick check of the panel
//   --bacdive_full  every BacDive strain with a NaCl maximum and a genome (BACDIVE_HARVEST), and the
//                   within-domain salt-in control's strains; in --nchunks tasks
// The random pool is built from the pilot strains in both cases, unless --random_panel gives one.
//
include { BACDIVE_HARVEST                          } from '../../../modules/local/bacdive/harvest/main'
include { BACDIVE_GENOMES as BACDIVE_GENOMES_PILOT } from '../../../modules/local/bacdive/genomes/main'
include { BACDIVE_GENOMES                          } from '../../../modules/local/bacdive/genomes/main'
include { BACDIVE_PROTEINS                         } from '../../../modules/local/bacdive/proteins/main'
include { RANDOM_POOL                              } from '../random_pool/main'

workflow BACDIVE_DATA {
    take:
    families        // families.yaml
    snapshot        // harvested BacDive batches to start from (bacdive_raw/), or []
    random_panel    // a given pool ([pool.hmm, pool_families.tsv]), or [] to build it

    main:
    pilot = file("${projectDir}/assets/bacdive/pilot_strains.tsv", checkIfExists: true)
    pilot_matrix = file("${projectDir}/assets/bacdive/feature_matrix.tsv", checkIfExists: true)
    n_pilot = pilot.readLines().drop(1).findAll { line -> line.trim() }.size()

    // the pilot's genomes: the analysis set by default, and the pool's screen
    if( !params.bacdive_full || !params.random_panel ) {
        BACDIVE_GENOMES_PILOT(channel.of(tuple([id: 'pilot'], pilot)))
        pilot_genomes = BACDIVE_GENOMES_PILOT.out.genomes.map { _meta, d -> d }.first()   // a value: every chunk reads it
    }

    // the analysis set
    if( params.bacdive_full ) {
        BACDIVE_HARVEST(snapshot)
        BACDIVE_GENOMES(BACDIVE_HARVEST.out.strains.map { f -> tuple([id: 'bacdive'], f) }
            .mix(BACDIVE_HARVEST.out.control_strains.map { f -> tuple([id: 'control'], f) }))
        // values: every chunk reads them
        genomes         = BACDIVE_GENOMES.out.genomes.filter { meta, _d -> meta.id == 'bacdive' }.map { _meta, d -> d }.first()
        control_genomes = BACDIVE_GENOMES.out.genomes.filter { meta, _d -> meta.id == 'control' }.map { _meta, d -> d }.first()
        strains         = BACDIVE_HARVEST.out.strains
        classes         = BACDIVE_HARVEST.out.classes
        control_strains = BACDIVE_HARVEST.out.control_strains
        nchunks         = params.nchunks
    } else {
        genomes         = pilot_genomes
        control_genomes = channel.value([])
        strains         = channel.value(pilot)
        classes         = channel.value([])
        control_strains = channel.value([])
        nchunks         = 1
    }

    // the random pool
    if( params.random_panel ) {
        pool_hmm      = random_panel.map { files -> files[0] }
        pool_families = random_panel.map { files -> files[1] }
    } else {
        RANDOM_POOL(families, pilot_genomes, pilot_matrix, n_pilot)
        pool_hmm      = RANDOM_POOL.out.hmm
        pool_families = RANDOM_POOL.out.families
    }

    // per genome: Prodigal proteins, size/GC, the pool search; each task every nchunks-th genome
    chunks = channel.of(0..(nchunks - 1)).map { i -> "${i}/${nchunks}".toString() }
    BACDIVE_PROTEINS(chunks, genomes, control_genomes, pool_hmm)

    emit:
    genomes         = genomes
    control_genomes = control_genomes
    pilot           = strains                                               // the strain table (pilot_strains.tsv)
    classes         = classes                                               // bacdive_classes.tsv, or []
    control_pilot   = control_strains                                       // the salt-in control's table, or []
    chunks          = BACDIVE_PROTEINS.out.chunk                            // [chunk, strains/, control_strains/]
    strains         = BACDIVE_PROTEINS.out.strains.flatten().toList()      // every strain's <bid>/
    control_strains = BACDIVE_PROTEINS.out.control_strains.flatten().toList()
    pool_families   = pool_families
}
