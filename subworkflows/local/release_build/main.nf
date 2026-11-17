//
// Build a release from the source release's snapshot: splits, alignments, HMMs, DIAMOND database,
// then the benchmark (calibration decides the cascade and report levels, the test reads score them).
// Every part is passed on as it is made; the steps that make a part of the release publish it into
// releases/<release> (conf/modules.config). RELEASE_DB is the one directory osmotool reads.
//
include { RELEASE_REFS       } from '../../../modules/local/release_refs/main'
include { SPLIT_REFS         } from '../../../modules/local/split_refs/main'
include { ALIGN_FAMILIES     } from '../../../modules/local/align_families/main'
include { BUILD_HMMS         } from '../../../modules/local/build_hmms/main'
include { BUILD_DIAMOND      } from '../../../modules/local/build_diamond/main'
include { RELEASE_DB as RELEASE_DB_CALIB } from '../../../modules/local/release/db/main'
include { RELEASE_DB                     } from '../../../modules/local/release/db/main'
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
    ALIGN_FAMILIES(SPLIT_REFS.out.refs, families, source)
    BUILD_HMMS(SPLIT_REFS.out.refs, ALIGN_FAMILIES.out.alignments, families)
    BUILD_DIAMOND(BUILD_HMMS.out.refs, families)
    refs          = BUILD_DIAMOND.out.refs
    release_files = BUILD_DIAMOND.out.release_files
    hmms          = BUILD_HMMS.out.hmms

    // calibration reads decide the cascade and report levels, on a database without either ...
    RELEASE_DB_CALIB(params.release, hmms, release_files, none)
    BENCHMARK_SIMULATE_CALIB(refs, families, 'calib')
    BENCHMARK_READS_CALIB(BENCHMARK_SIMULATE_CALIB.out.reads, RELEASE_DB_CALIB.out.db, 'score')
    BENCHMARK_METRICS_CALIB(BENCHMARK_SIMULATE_CALIB.out.reads, BENCHMARK_READS_CALIB.out.results, none, families, 'calib')
    cascade = BENCHMARK_METRICS_CALIB.out.cascade

    // ... the release's database has them, and the held-out test reads score it
    RELEASE_DB(params.release, hmms, release_files, cascade)
    db = RELEASE_DB.out.db
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
    db            = db                            // the release as osmotool reads it
    refs          = refs
    alignments    = BUILD_HMMS.out.alignments
    hmms          = hmms
    release_files = release_files                 // the DIAMOND database and its files
    cascade       = cascade                       // cascade config and report levels
    done          = BENCHMARK_METRICS_TEST.out.metrics.mix(BENCHMARK_ANNOTATESCORE.out.annotate).collect()
}
