//
// After the gate: BacDive annotation with the frozen release (nchunks tasks, each on the proteins of the
// same chunk), the strain x family matrix, the phenotype validation, and the NaCl model.
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

    main:
    BACDIVE_ANNOTATE(chunks, genomes, control_genomes, db)
    annotate         = BACDIVE_ANNOTATE.out.annotate.flatten().toList()
    control_annotate = BACDIVE_ANNOTATE.out.control_annotate.flatten().toList()

    BACDIVE_MATRIX(genomes, pilot, annotate)
    matrix = BACDIVE_MATRIX.out.matrix
    BACDIVE_ANALYZE(matrix, pool_families, strains, classes, control_pilot, control_strains, control_annotate)
    NACL_MODEL(db, matrix, annotate, strains, BACDIVE_ANALYZE.out.family_ranking)

    emit:
    done = BACDIVE_ANALYZE.out.reports.mix(NACL_MODEL.out.model).collect()
}
