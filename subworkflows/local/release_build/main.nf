//
// Build a release from the source release's snapshot: splits, alignments, HMMs, DIAMOND database,
// then the benchmark (calibration decides the cascade and report levels, the test reads score them).
// Every part is passed on as it is made; the steps that make a part of the release publish it into
// releases/<release> (conf/modules.config).
//
include { family_metas; alignment_plan; sha256 } from '../utils_osmo_refdb_pipeline/main'
include { RELEASE_REFS       } from '../../../modules/local/release_refs/main'
include { SPLIT_REFS         } from '../../../modules/local/split_refs/main'
include { MAFFT_ALIGN        } from '../../../modules/local/mafft/align/main'
include { TRIMAL             } from '../../../modules/local/trimal/main'
include { BUILD_HMMS         } from '../../../modules/local/build_hmms/main'
include { BUILD_DIAMOND      } from '../../../modules/local/build_diamond/main'
include { BENCHMARK_SIMULATE as BENCHMARK_SIMULATE_CALIB } from '../../../modules/local/benchmark/simulate/main'
include { BENCHMARK_SIMULATE as BENCHMARK_SIMULATE_TEST  } from '../../../modules/local/benchmark/simulate/main'
include { BENCHMARK_READS    as BENCHMARK_READS_CALIB    } from '../../../modules/local/benchmark/reads/main'
include { BENCHMARK_READS    as BENCHMARK_READS_TEST     } from '../../../modules/local/benchmark/reads/main'
include { BENCHMARK_METRICS  as BENCHMARK_METRICS_CALIB  } from '../../../modules/local/benchmark/metrics/main'
include { BENCHMARK_METRICS  as BENCHMARK_METRICS_TEST   } from '../../../modules/local/benchmark/metrics/main'
include { BENCHMARK_ANNOTATEPREP  } from '../../../modules/local/benchmark/annotateprep/main'
include { OSMOTOOL_ANNOTATE as OSMOTOOL_ANNOTATE_BENCHMARK } from '../../../modules/local/osmotool/annotate/main'
include { BENCHMARK_ANNOTATESCORE } from '../../../modules/local/benchmark/annotatescore/main'

workflow RELEASE_BUILD {
    take:
    families    // families.yaml
    source      // the source release's snapshot (--reuse_refs_from)

    main:
    none = channel.value([])

    RELEASE_REFS(families, source)
    SPLIT_REFS(RELEASE_REFS.out.refs, families)

    // ---- alignments, one family at a time: the families in families.yaml that have a train set
    metas = families.flatMap { f -> family_metas(f) }.map { meta -> tuple(meta.id, meta) }
    train = SPLIT_REFS.out.train.flatten().map { f -> tuple(f.name - '.positive.train.faa', f) }
    planned = metas.join(train)
        .combine(source)
        .map { _id, meta, f, src -> tuple(meta, f, alignment_plan(meta, f, src), src) }
    plan = planned.branch { _meta, _f, how, _src ->
        reuse: how == 'reuse'
        align: how == 'align'
        skip:  true
    }
    MAFFT_ALIGN(plan.align.map { meta, f, _how, _src -> tuple(meta, f) })
    TRIMAL(MAFFT_ALIGN.out.alignment)
    reused = plan.reuse.flatMap { meta, _f, _how, src ->
        ['aln.fasta', 'trimmed.fasta'].collect { ext -> src.resolve("alignments/${meta.id}.${ext}") }
    }
    // what each family's alignment is, for the next build to check (see alignment_plan)
    manifest = planned
        .map { meta, f, how, src ->
            def what = how == 'reuse' ? "reused from ${src.name}" : how == 'align' ? 'aligned' : 'none (fewer than 2 sequences)'
            "${meta.id}\t${sha256(f)}\t${meta.trim_gt}\t${what}"
        }
        .collectFile(name: 'manifest.tsv', newLine: true, sort: true, seed: "family\ttrain_sha256\ttrim_gt\talignment")
    // the source's HMM-decoy alignments: BUILD_HMMS (05b) reuses them when the decoy's input is unchanged
    decoy_alignments = source.flatMap { src -> files("${src}/alignments/*.aln") }
    alignments = MAFFT_ALIGN.out.alignment.mix(TRIMAL.out.trimmed).map { _meta, f -> f }
        .mix(reused, decoy_alignments, manifest)
        .collect()

    BUILD_HMMS(SPLIT_REFS.out.refs, alignments, families)
    BUILD_DIAMOND(BUILD_HMMS.out.refs, families)
    refs          = BUILD_DIAMOND.out.refs
    release_files = BUILD_DIAMOND.out.release_files
    hmms          = BUILD_HMMS.out.hmms

    // ---- the benchmark. osmotool reads one directory (DATABASE): each osmotool step is given the
    // DIAMOND files, the cascade and hmms/ as one list, staged as that directory
    // calibration reads decide the cascade and report levels, on a database without either ...
    db_calib = release_files.mix(hmms).collect()
    BENCHMARK_SIMULATE_CALIB(refs, families, 'calib')
    BENCHMARK_READS_CALIB(BENCHMARK_SIMULATE_CALIB.out.reads, db_calib, 'score')
    BENCHMARK_METRICS_CALIB(BENCHMARK_SIMULATE_CALIB.out.reads, BENCHMARK_READS_CALIB.out.results, none, families, 'calib')
    cascade = BENCHMARK_METRICS_CALIB.out.cascade

    // ... the release's database has them, and the held-out test reads score it
    db = release_files.mix(cascade, hmms).collect()
    BENCHMARK_SIMULATE_TEST(refs, families, 'test')
    BENCHMARK_READS_TEST(BENCHMARK_SIMULATE_TEST.out.reads, db, 'on')
    BENCHMARK_METRICS_TEST(BENCHMARK_SIMULATE_TEST.out.reads, BENCHMARK_READS_TEST.out.results, cascade, families, 'test')

    // osmotool annotate on held-out full-length proteins, once per method (ext.args: --method)
    BENCHMARK_ANNOTATEPREP(refs, families)
    eval_set = channel.of('hmm', 'diamond')
        .combine(BENCHMARK_ANNOTATEPREP.out.eval_set)
        .map { method, dummy, proteins -> tuple([id: "annotate_${method}", method: method], dummy, proteins) }
    OSMOTOOL_ANNOTATE_BENCHMARK(eval_set, db)
    BENCHMARK_ANNOTATESCORE(refs, families,
        OSMOTOOL_ANNOTATE_BENCHMARK.out.results.map { _meta, files -> files }.flatten().collect())

    emit:
    db            = db                            // the release as osmotool reads it: [files..., hmms/]
    refs          = refs
    alignments    = BUILD_HMMS.out.alignments
    hmms          = hmms
    release_files = release_files                 // the DIAMOND database and its files
    cascade       = cascade                       // cascade config and report levels
    done          = BENCHMARK_METRICS_TEST.out.metrics.mix(BENCHMARK_ANNOTATESCORE.out.annotate).collect()
}
