//
// The release's alignments and HMMs, one family (or decoy) per task:
//   families  train set -> MAFFT_ALIGN -> TRIMAL -> HMMER_HMMBUILD   (or PFAM_FETCH, for a family with a
//             curated pfam_model); fewer than 2 train sequences: no alignment, no model
//   decoys    HMM_DECOYSELECT -> CDHIT -> MAFFT_ALIGN -> TRIMAL -> HMMER_HMMBUILD
//   scoring   HMMER_HMMSEARCH: each model on its family's calibration sets, each decoy-aware decoy on
//             its competitors'
//   then      HMM_CALIBRATE: GA cutoffs, the QC scorecard, one pressed database
//
include { family_metas; decoy_metas; scoring_sets; refs_kind; alignment_plan } from '../utils_osmo_refdb_pipeline/main'
include { MAFFT_ALIGN                              } from '../../../modules/local/mafft/align/main'
include { MAFFT_ALIGN    as MAFFT_ALIGN_DECOY      } from '../../../modules/local/mafft/align/main'
include { TRIMAL                                   } from '../../../modules/local/trimal/main'
include { TRIMAL         as TRIMAL_DECOY           } from '../../../modules/local/trimal/main'
include { HMMER_HMMBUILD                           } from '../../../modules/local/hmmer/hmmbuild/main'
include { HMMER_HMMBUILD as HMMER_HMMBUILD_DECOY   } from '../../../modules/local/hmmer/hmmbuild/main'
include { PFAM_FETCH                               } from '../../../modules/local/pfam/fetch/main'
include { HMM_DECOYSELECT                          } from '../../../modules/local/hmm/decoyselect/main'
include { CDHIT          as CDHIT_DECOY            } from '../../../modules/local/cdhit/main'
include { HMMER_HMMSEARCH                          } from '../../../modules/local/hmmer/hmmsearch/main'
include { HMM_CALIBRATE                            } from '../../../modules/local/hmm/calibrate/main'

workflow HMM_BUILD {
    take:
    families    // families.yaml
    refs        // SPLIT_REFS: the refs directory
    train       // SPLIT_REFS: every family's positive.train set
    sets        // SPLIT_REFS: every refs/*.faa

    main:
    metas = families.flatMap { f -> family_metas(f) }

    // ---- families: alignment (see alignment_plan)
    planned = metas.map { meta -> tuple(meta.id, meta) }
        .join(train.flatten().map { f -> tuple(f.name - '.positive.train.faa', f) })
        .map { _id, meta, f -> tuple(meta, f, alignment_plan(f)) }
    MAFFT_ALIGN(planned.filter { _meta, _f, how -> how == 'align' }.map { meta, f, _how -> tuple(meta, f) })
    TRIMAL(MAFFT_ALIGN.out.alignment)
    trimmed = TRIMAL.out.trimmed

    // ---- families: the model, built from the trimmed alignment or Pfam's own
    HMMER_HMMBUILD(trimmed.filter { meta, _f -> !meta.pfam_model })
    PFAM_FETCH(metas.filter { meta -> meta.pfam_model })
    models = HMMER_HMMBUILD.out.hmm.mix(PFAM_FETCH.out.hmm)

    // ---- decoys: input set, then clustering, alignment and trimming
    HMM_DECOYSELECT(refs, families)
    decoys = families.flatMap { f -> decoy_metas(f) }.map { meta -> tuple(meta.id, meta) }
        .join(HMM_DECOYSELECT.out.inputs.flatten().map { f -> tuple(f.name - '.hmm_decoy.faa', f) })
        .map { _id, meta, f -> tuple(meta, f) }
    CDHIT_DECOY(decoys)
    MAFFT_ALIGN_DECOY(CDHIT_DECOY.out.clustered)
    TRIMAL_DECOY(MAFFT_ALIGN_DECOY.out.alignment)
    HMMER_HMMBUILD_DECOY(TRIMAL_DECOY.out.trimmed)
    decoy_models = HMMER_HMMBUILD_DECOY.out.hmm

    // ---- scoring: each model on its family's calibration sets ...
    by_family = sets.flatten()
        .filter { f -> f.size() > 0 }
        .map { f -> refs_kind(f) }
        .groupTuple()
        .map { fam, kinds, files -> tuple(fam, [kinds, files].transpose().collectEntries()) }
    family_jobs = models.map { meta, hmm -> tuple(meta.id, hmm) }
        .join(by_family)
        .flatMap { fam, hmm, kinds ->
            def chosen = scoring_sets(fam, kinds)
            ['positive', 'negative'].findAll { kind -> chosen[kind] }
                .collect { kind -> tuple([id: "${fam}.${kind}"], hmm, chosen[kind]) }
        }
    // ... and each decoy-aware decoy on its competitors' positive.calib and negative.train
    decoy_jobs = decoy_models.filter { meta, _hmm -> meta.decoy_aware }
        .combine(by_family.toList().map { list -> [list.collectEntries()] })
        .flatMap { meta, hmm, all ->
            meta.competes_with.collectMany { fam ->
                [positive: 'positive.calib', negative: 'negative.train']
                    .findAll { _kind, split -> all[fam]?.containsKey(split) }
                    .collect { kind, split -> tuple([id: "${fam}.${meta.id}.${kind}"], hmm, all[fam][split]) }
            }
        }
    HMMER_HMMSEARCH(family_jobs.mix(decoy_jobs))

    // ---- GA cutoffs, scorecard, the pressed database
    HMM_CALIBRATE(
        models.mix(decoy_models).map { _meta, hmm -> hmm }.collect(),
        HMMER_HMMSEARCH.out.tblout.map { _meta, f -> f }.collect(),
        refs, families)

    emit:
    hmms       = HMM_CALIBRATE.out.hmms                  // the release's hmms/
    qc         = HMM_CALIBRATE.out.qc                    // qc_scorecard.tsv
    // the release's alignments/ (each step publishes its own, conf/modules.config)
    alignments = MAFFT_ALIGN.out.alignment.mix(TRIMAL.out.trimmed, MAFFT_ALIGN_DECOY.out.alignment, TRIMAL_DECOY.out.trimmed)
                     .map { _meta, f -> f }.collect()
    // the decoys' input and clustered sets, for the release's refs/
    decoy_refs = HMM_DECOYSELECT.out.inputs.flatten()
                     .mix(CDHIT_DECOY.out.clustered.mix(CDHIT_DECOY.out.clusters).map { _meta, f -> f })
                     .collect().ifEmpty([])
}
