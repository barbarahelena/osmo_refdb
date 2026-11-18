//
// The random-panel pool: the non-osmotic Pfam families whose random panels are the BacDive analysis's
// null baseline. Part of the analysis's design, so it is built from fixed inputs and does not depend on
// the release: Pfam-A 39.0, the 43 pilot strains' proteomes as the screen, and the pilot's osmo feature
// matrix (assets/bacdive/) to match the panel's prevalence; the seeds are fixed (conf/modules.config).
//   POOL_PFAM -> POOL_CANDIDATES -> hmmsearch --cut_ga on POOL_PROTEOMES -> POOL_BUILD
//
include { POOL_PFAM                         } from '../../../modules/local/pool/pfam/main'
include { POOL_PROTEOMES                    } from '../../../modules/local/pool/proteomes/main'
include { POOL_CANDIDATES                   } from '../../../modules/local/pool/candidates/main'
include { HMMER_HMMSEARCH as POOL_SCREEN    } from '../../../modules/local/hmmer/hmmsearch/main'
include { POOL_BUILD                        } from '../../../modules/local/pool/build/main'

workflow RANDOM_POOL {
    take:
    families        // families.yaml: its Pfam accessions are never pool families
    genomes         // the pilot strains' genomes (BACDIVE_GENOMES)
    osmo_matrix     // the pilot's osmo feature matrix
    n_genomes       // how many screening genomes

    main:
    POOL_PFAM('39.0')
    POOL_PROTEOMES(genomes)
    POOL_CANDIDATES(POOL_PFAM.out.pfam, families)
    POOL_SCREEN(POOL_CANDIDATES.out.hmm.combine(POOL_PROTEOMES.out.proteomes)
        .map { hmm, proteomes -> tuple([id: 'screen'], hmm, proteomes) })
    POOL_BUILD(POOL_PFAM.out.pfam, POOL_CANDIDATES.out.table,
               POOL_SCREEN.out.tblout.map { _meta, f -> f }, osmo_matrix, n_genomes)

    emit:
    hmm      = POOL_BUILD.out.hmm.first()         // pool.hmm (a value: every chunk searches it)
    families = POOL_BUILD.out.families.first()    // pool_families.tsv
}
