//
// The release's alignments and HMMs, one family (or decoy) per task:
//   families  train set -> MAFFT_ALIGN -> TRIMAL -> HMMER_HMMBUILD   (or the source's alignment, reused;
//             or PFAM_FETCH, for a family with a curated pfam_model)
//   decoys    HMM_DECOYSELECT -> CDHIT -> MAFFT_ALIGN -> TRIMAL -> HMMER_HMMBUILD   (or reused)
//   scoring   HMMER_HMMSEARCH: each model on its family's calibration sets, each decoy-aware decoy on
//             its competitors'
//   then      HMM_CALIBRATE: GA cutoffs, the QC scorecard, one pressed database
// What can be reused from the source release is decided per family in the workflow
// (alignment_plan, decoy_plan), and recorded in alignments/manifest.tsv.
//
include { family_metas; decoy_metas; scoring_sets; refs_kind; alignment_plan; decoy_plan; source_decoy_files; sha256 } from '../utils_osmo_refdb_pipeline/main'
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
include { ALIGNMENTS_RECORD                        } from '../../../modules/local/alignments/record/main'

workflow HMM_BUILD {
    take:
    families    // families.yaml
    source      // the source release's snapshot
    refs        // SPLIT_REFS: the refs directory
    train       // SPLIT_REFS: every family's positive.train set
    sets        // SPLIT_REFS: every refs/*.faa

    main:
    metas = families.flatMap { f -> family_metas(f) }

    // ---- families: alignment (see alignment_plan)
    planned = metas.map { meta -> tuple(meta.id, meta) }
        .join(train.flatten().map { f -> tuple(f.name - '.positive.train.faa', f) })
        .combine(source)
        .map { _id, meta, f, src -> tuple(meta, f, alignment_plan(meta, f, src), src) }
    plan = planned.branch { _meta, _f, how, _src ->
        reuse: how == 'reuse'
        align: how == 'align'
        skip:  true
    }
    MAFFT_ALIGN(plan.align.map { meta, f, _how, _src -> tuple(meta, f) })
    TRIMAL(MAFFT_ALIGN.out.alignment)
    reused = plan.reuse.map { meta, _f, _how, src ->
        tuple(meta, src.resolve("alignments/${meta.id}.aln.fasta"), src.resolve("alignments/${meta.id}.trimmed.fasta"))
    }
    trimmed = TRIMAL.out.trimmed.mix(reused.map { meta, _raw, trim -> tuple(meta, trim) })

    // ---- families: the model, built from the trimmed alignment or Pfam's own
    HMMER_HMMBUILD(trimmed.filter { meta, _f -> !meta.pfam_model })
    PFAM_FETCH(metas.filter { meta -> meta.pfam_model })
    models = HMMER_HMMBUILD.out.hmm.mix(PFAM_FETCH.out.hmm)

    // ---- decoys: input set, then clustering, alignment and trimming unless reused (see decoy_plan)
    HMM_DECOYSELECT(refs, families)
    decoys = families.flatMap { f -> decoy_metas(f) }.map { meta -> tuple(meta.id, meta) }
        .join(HMM_DECOYSELECT.out.inputs.flatten().map { f -> tuple(f.name - '.hmm_decoy.faa', f) })
        .combine(source)
        .map { _id, meta, f, src -> tuple(meta, f, decoy_plan(meta, f, src), src) }
    dplan = decoys.branch { _meta, _f, how, _src ->
        reuse: how == 'reuse'
        align: true
    }
    CDHIT_DECOY(dplan.align.map { meta, f, _how, _src -> tuple(meta, f) })
    MAFFT_ALIGN_DECOY(CDHIT_DECOY.out.clustered)
    TRIMAL_DECOY(MAFFT_ALIGN_DECOY.out.alignment)
    dreused = dplan.reuse.map { meta, _f, _how, src -> tuple(meta, source_decoy_files(src, meta.id)) }
    HMMER_HMMBUILD_DECOY(TRIMAL_DECOY.out.trimmed.mix(dreused.map { meta, old -> tuple(meta, old.trimmed) }))
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

    // ---- the alignments and decoy sets as the release keeps them, and what each one is
    manifest = planned.map { meta, f, how, src -> tuple(meta, f, how, src.name) }
        .mix(decoys.map { meta, f, how, src -> tuple(meta, f, how, src.name) })
        .map { meta, f, how, src_name ->
            def what = how == 'reuse' ? "reused from ${src_name}" : how == 'align' ? 'aligned' : 'none (fewer than 2 sequences)'
            "${meta.id}\t${sha256(f)}\t${meta.trim_gt}\t${what}"
        }
        .collectFile(name: 'manifest.tsv', newLine: true, sort: true, seed: "family\tinput_sha256\ttrim_gt\talignment")
    alignment_files = MAFFT_ALIGN.out.alignment.mix(TRIMAL.out.trimmed, MAFFT_ALIGN_DECOY.out.alignment, TRIMAL_DECOY.out.trimmed)
        .map { _meta, f -> f }
        .mix(reused.flatMap { _meta, raw, trim -> [raw, trim] },
             dreused.flatMap { _meta, old -> [old.raw, old.trimmed] },
             manifest)
        .collect()
    decoy_refs = HMM_DECOYSELECT.out.inputs.flatten()
        .mix(CDHIT_DECOY.out.clustered.mix(CDHIT_DECOY.out.clusters).map { _meta, f -> f },
             dreused.flatMap { _meta, old -> [old.clustered, old.clusters].findAll { f -> f.exists() } })
        .collect()
        .ifEmpty([])
    ALIGNMENTS_RECORD(alignment_files, decoy_refs)

    emit:
    hmms       = HMM_CALIBRATE.out.hmms                  // the release's hmms/
    qc         = HMM_CALIBRATE.out.qc                    // qc_scorecard.tsv
    alignments = ALIGNMENTS_RECORD.out.alignments        // the release's alignments/
    decoy_refs = ALIGNMENTS_RECORD.out.decoy_refs        // the decoys' input and clustered sets
}
