//
// The release's DIAMOND database and its files:
//   DIAMOND_EXTRAREFS   decoy and fused-ORF references
//   DIAMOND_TRAINREFS   the sequences the database is built from (<name>.train_refs.faa)
//   DIAMOND_MAKEDB      <name>.dmnd
//   DIAMOND_BLASTP      per family, its calibration sets against the database
//   DIAMOND_CALIBRATE   the per-family bit-score cutoffs and the excluded-family lists
//
include { diamond_sets; refs_kind } from '../utils_osmo_refdb_pipeline/main'
include { DIAMOND_EXTRAREFS } from '../../../modules/local/diamond/extrarefs/main'
include { DIAMOND_TRAINREFS } from '../../../modules/local/diamond/trainrefs/main'
include { DIAMOND_MAKEDB    } from '../../../modules/local/diamond/makedb/main'
include { DIAMOND_BLASTP    } from '../../../modules/local/diamond/blastp/main'
include { DIAMOND_CALIBRATE } from '../../../modules/local/diamond/calibrate/main'

workflow DIAMOND_BUILD {
    take:
    families    // families.yaml
    refs        // SPLIT_REFS: the refs directory
    train       // SPLIT_REFS: every family's positive.train set
    sets        // SPLIT_REFS: every refs/*.faa

    main:
    name = 'osmo_refdb'   // the fixed file names osmotool looks for

    DIAMOND_EXTRAREFS(refs, families)
    extra = DIAMOND_EXTRAREFS.out.refs.flatten().collect().ifEmpty([])
    DIAMOND_TRAINREFS(name, train, extra)
    DIAMOND_MAKEDB(DIAMOND_TRAINREFS.out.fasta.map { f -> tuple([id: name], f) })
    db = DIAMOND_MAKEDB.out.db.map { _meta, dmnd -> dmnd }

    // per family, its calibration sets (diamond_sets) against the whole database
    jobs = sets.flatten()
        .map { f -> refs_kind(f) }
        .groupTuple()
        .flatMap { fam, kinds, files ->
            diamond_sets([kinds, files].transpose().collectEntries())
                .collect { kind, f -> tuple([id: "${fam}.${kind}"], f) }
        }
    DIAMOND_BLASTP(jobs, db)
    DIAMOND_CALIBRATE(name, refs, DIAMOND_BLASTP.out.hits.map { _meta, f -> f }.collect().ifEmpty([]), families)

    emit:
    extra         = extra                                            // the decoy and fused-ORF references
    release_files = DIAMOND_TRAINREFS.out.fasta.mix(db, DIAMOND_CALIBRATE.out.cutoffs, DIAMOND_CALIBRATE.out.excluded)
                        .flatten().collect()                         // the release's DIAMOND files
}
