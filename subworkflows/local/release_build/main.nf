//
// Build a release from the source release's snapshot: splits, alignments, HMMs, DIAMOND database,
// then the benchmark (calibration decides the cascade and report levels, the test reads score them).
//
include { under_repo         } from '../utils_osmo_refdb_pipeline/main'
include { RELEASE_REFS       } from '../../../modules/local/release_refs/main'
include { SPLIT_REFS         } from '../../../modules/local/split_refs/main'
include { ALIGN_FAMILIES     } from '../../../modules/local/align_families/main'
include { BUILD_HMMS         } from '../../../modules/local/build_hmms/main'
include { BUILD_DIAMOND      } from '../../../modules/local/build_diamond/main'
include { SIMULATE_READS as SIMULATE_CALIB } from '../../../modules/local/simulate_reads/main'
include { SIMULATE_READS as SIMULATE_TEST  } from '../../../modules/local/simulate_reads/main'
include { RUN_BENCHMARK as RUN_CALIB       } from '../../../modules/local/run_benchmark/main'
include { RUN_BENCHMARK as RUN_TEST        } from '../../../modules/local/run_benchmark/main'
include { COMPUTE_METRICS as METRICS_CALIB } from '../../../modules/local/compute_metrics/main'
include { COMPUTE_METRICS as METRICS_TEST  } from '../../../modules/local/compute_metrics/main'
include { BENCHMARK_ANNOTATE } from '../../../modules/local/benchmark_annotate/main'
include { ASSEMBLE_RELEASE   } from '../../../modules/local/assemble_release/main'

workflow RELEASE_BUILD {
    main:
    families = channel.value(file(params.families_file))
    source   = channel.value(file(under_repo(params.reuse_refs_from), checkIfExists: true))
    none     = channel.value([])

    RELEASE_REFS(families, source)
    SPLIT_REFS(RELEASE_REFS.out.refs, families)
    ALIGN_FAMILIES(SPLIT_REFS.out.refs, families, source)
    BUILD_HMMS(SPLIT_REFS.out.refs, ALIGN_FAMILIES.out.alignments, families)
    BUILD_DIAMOND(BUILD_HMMS.out.refs, families)
    refs = BUILD_DIAMOND.out.refs
    rel  = BUILD_DIAMOND.out.release_files
    hmms = BUILD_HMMS.out.hmms

    // calibration reads decide the cascade and report levels ...
    SIMULATE_CALIB(refs, families, 'calib')
    RUN_CALIB(SIMULATE_CALIB.out.reads, hmms, rel, none, 'score')
    METRICS_CALIB(SIMULATE_CALIB.out.reads, RUN_CALIB.out.results, none, families, 'calib')

    // ... the held-out test reads score them
    SIMULATE_TEST(refs, families, 'test')
    RUN_TEST(SIMULATE_TEST.out.reads, hmms, rel, METRICS_CALIB.out.cascade, 'on')
    METRICS_TEST(SIMULATE_TEST.out.reads, RUN_TEST.out.results, METRICS_CALIB.out.cascade, families, 'test')
    BENCHMARK_ANNOTATE(refs, hmms, rel, METRICS_CALIB.out.cascade, families)

    ASSEMBLE_RELEASE(refs, BUILD_HMMS.out.alignments, hmms, rel, BUILD_HMMS.out.qc,
                     METRICS_CALIB.out.cascade, METRICS_CALIB.out.metrics,
                     METRICS_TEST.out.metrics, BENCHMARK_ANNOTATE.out.annotate)

    emit:
    release = ASSEMBLE_RELEASE.out
}
