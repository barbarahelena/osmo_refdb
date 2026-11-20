//
// The BacDive strains, their genomes, and per genome its Prodigal proteins and random-pool search.
// Release-independent: runs alongside the build.
//   default         the 43 pilot strains (assets/bacdive/pilot_strains.tsv, or --bacdive_pilot): a quick check of the panel
//   --bacdive_full  every BacDive strain with a NaCl maximum and a genome (BACDIVE_HARVEST), and the
//                   within-domain salt-in control's strains; in --nchunks tasks
// The random-panel pool is the same in both cases.
//
include { BACDIVE_HARVEST                          } from '../../../modules/local/bacdive/harvest/main'
include { BACDIVE_GENOMES as BACDIVE_GENOMES_PILOT } from '../../../modules/local/bacdive/genomes/main'
include { BACDIVE_GENOMES                          } from '../../../modules/local/bacdive/genomes/main'
include { BACDIVE_PROTEINS                         } from '../../../modules/local/bacdive/proteins/main'
include { POOL_PFAM                                } from '../../../modules/local/pool/pfam/main'
include { POOL_EXTRACT                             } from '../../../modules/local/pool/extract/main'

workflow BACDIVE_DATA {
    take:
    snapshot        // harvested BacDive batches to start from (bacdive_raw/), or []

    main:
    pilot = file(params.bacdive_pilot ?: "${projectDir}/assets/bacdive/pilot_strains.tsv", checkIfExists: true)

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
        BACDIVE_GENOMES_PILOT(channel.of(tuple([id: 'pilot'], pilot)))
        genomes         = BACDIVE_GENOMES_PILOT.out.genomes.map { _meta, d -> d }.first()   // a value: every chunk reads it
        control_genomes = channel.value([])
        strains         = channel.value(pilot)
        classes         = channel.value([])
        control_strains = channel.value([])
        nchunks         = 1
    }

    // the random-panel pool, the analysis's null baseline: fixed, the 174 non-osmotic Pfam families of
    // assets/random_panel/pool_families.tsv (drawn from Pfam-A 39.0, screened on the 43 pilot proteomes,
    // stratified to the osmo panel's prevalence), their models taken from Pfam-A 39.0
    pool_families = channel.value(file("${projectDir}/assets/random_panel/pool_families.tsv", checkIfExists: true))
    POOL_PFAM('39.0')
    POOL_EXTRACT(POOL_PFAM.out.pfam, pool_families)
    pool_hmm = POOL_EXTRACT.out.hmm   // a value (its inputs are): every chunk searches it

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
