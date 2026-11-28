//
// On the finished release: BacDive annotation with the frozen release (each task on the proteins of the same
// chunk), the strain x family matrix, the phenotype validation, and, with --bacdive_full, the NaCl
// model (41 pilot strains are too few to train it on).
//
include { BACDIVE_ANNOTATE } from '../../../modules/local/bacdive/annotate/main'
include { BACDIVE_MATRIX   } from '../../../modules/local/bacdive/matrix/main'
include { BACDIVE_ANALYZE  } from '../../../modules/local/bacdive/analyze/main'
include { NACL_MODEL       } from '../../../modules/local/nacl/model/main'

workflow BACDIVE_ANALYSIS {
    take:
    db                // the release as osmotool reads it
    chunks            // BACDIVE_DATA: [chunk, strains/, control_strains/]
    genomes           //               the genomes
    control_genomes   //               the salt-in control's genomes
    pilot             //               pilot_strains.tsv
    classes           //               bacdive_classes.tsv
    control_pilot     //               the control's pilot_strains.tsv
    strains           //               every strain's <bid>/ (pool search, proteins)
    control_strains   //               the control's
    pool_families     // the random-panel pool, pool_families.tsv
    annotate_levels   // osmo_refdb.annotate_report_levels.tsv (the NaCl model's features)

    main:
    BACDIVE_ANNOTATE(chunks, genomes, control_genomes, db)
    annotate         = BACDIVE_ANNOTATE.out.annotate.flatten().toList()
    control_annotate = BACDIVE_ANNOTATE.out.control_annotate.flatten().toList()

    BACDIVE_MATRIX(genomes, pilot, annotate)
    matrix = BACDIVE_MATRIX.out.matrix
    BACDIVE_ANALYZE(matrix, pool_families, strains, classes, control_pilot, control_strains, control_annotate)
    done = BACDIVE_ANALYZE.out.reports
    if( params.bacdive_full ) {
        // the osmotool version the training genomes were annotated with, recorded in the model's schema
        osmotool_version = BACDIVE_ANNOTATE.out.versions_osmotool.map { _proc, _tool, v -> v }.unique().toList()
            .map { vs -> vs.size() == 1 ? vs[0] : error("BACDIVE_ANNOTATE ran more than one osmotool version: ${vs}") }
        NACL_MODEL(params.release, osmotool_version, db, matrix, annotate, strains, BACDIVE_ANALYZE.out.family_ranking,
            annotate_levels)
        done = done.mix(NACL_MODEL.out.model)
    }

    emit:
    done = done.collect()
}
