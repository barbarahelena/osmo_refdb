//
// Build a release from scratch: references (UniProt, RefSeq), splits, alignments, HMMs, DIAMOND database,
// then the benchmark (calibration decides the cascade and report levels, the test reads score them).
// Every part is passed on as it is made; the steps that make a part of the release publish it into
// <outdir> (conf/modules.config).
//
include { HMM_BUILD          } from '../hmm_build/main'
include { REFS_BUILD         } from '../refs_build/main'
include { SPLIT_REFS         } from '../../../modules/local/split_refs/main'
include { DIAMOND_BUILD      } from '../diamond_build/main'
include { BENCHMARK_SIMULATE as BENCHMARK_SIMULATE_CALIB } from '../../../modules/local/benchmark/simulate/main'
include { BENCHMARK_SIMULATE as BENCHMARK_SIMULATE_TEST  } from '../../../modules/local/benchmark/simulate/main'
include { BENCHMARK_READS    as BENCHMARK_READS_CALIB    } from '../../../modules/local/benchmark/reads/main'
include { BENCHMARK_READS    as BENCHMARK_READS_TEST     } from '../../../modules/local/benchmark/reads/main'
include { BENCHMARK_METRICS  as BENCHMARK_METRICS_CALIB  } from '../../../modules/local/benchmark/metrics/main'
include { BENCHMARK_METRICS  as BENCHMARK_METRICS_TEST   } from '../../../modules/local/benchmark/metrics/main'
include { BENCHMARK_ANNOTATEPREP  } from '../../../modules/local/benchmark/annotateprep/main'
include { OSMOTOOL_ANNOTATE as OSMOTOOL_ANNOTATE_BENCHMARK } from '../../../modules/local/osmotool/annotate/main'
include { BENCHMARK_ANNOTATESCORE } from '../../../modules/local/benchmark/annotatescore/main'
include { BENCHMARK_CATALOGGENES } from '../../../modules/local/benchmark/cataloggenes/main'
include { OSMOTOOL_BUILDCATALOG as OSMOTOOL_BUILDCATALOG_BENCHMARK } from '../../../modules/local/osmotool/buildcatalog/main'
include { BENCHMARK_CATALOGMAP   } from '../../../modules/local/benchmark/catalogmap/main'
include { BENCHMARK_CATALOGSCORE } from '../../../modules/local/benchmark/catalogscore/main'

workflow RELEASE_BUILD {
    take:
    families    // families.yaml

    main:
    none = channel.value([])

    REFS_BUILD(families, true)
    SPLIT_REFS(REFS_BUILD.out.refs, families)

    // alignments and HMMs, one family (or decoy) at a time
    HMM_BUILD(families, SPLIT_REFS.out.refs, SPLIT_REFS.out.train, SPLIT_REFS.out.faa)
    hmms = HMM_BUILD.out.hmms

    // the DIAMOND database and its files
    DIAMOND_BUILD(families, SPLIT_REFS.out.refs, SPLIT_REFS.out.train, SPLIT_REFS.out.faa)
    release_files = DIAMOND_BUILD.out.release_files

    // the release's refs/: the split sets, the DIAMOND decoy and fused-ORF references, the HMM decoys'
    // sets; given to the benchmark and provenance as one list, staged as refs/
    refs = SPLIT_REFS.out.files.mix(DIAMOND_BUILD.out.extra, HMM_BUILD.out.decoy_refs).flatten().collect()

    // ---- the benchmark. osmotool reads one directory (DATABASE): each osmotool step is given the
    // DIAMOND files, the cascade and hmms/ as one list, staged as that directory
    // calibration reads decide the cascade and report levels, on a database without either ...
    db_calib = release_files.mix(hmms).collect()
    BENCHMARK_SIMULATE_CALIB(refs, families, 'calib', none)
    BENCHMARK_READS_CALIB(BENCHMARK_SIMULATE_CALIB.out.reads, db_calib, 'score')
    BENCHMARK_METRICS_CALIB(BENCHMARK_SIMULATE_CALIB.out.reads, BENCHMARK_READS_CALIB.out.results, none, families, 'calib')
    cascade = BENCHMARK_METRICS_CALIB.out.cascade

    // ... the release's database has them, and the held-out test reads score it
    db = release_files.mix(cascade, hmms).collect()
    // a real metagenome as extra background (--background_r1/_r2), test reads only
    background = params.background_r1 ? channel.value([file(params.background_r1, checkIfExists: true), file(params.background_r2, checkIfExists: true)]) : none
    BENCHMARK_SIMULATE_TEST(refs, families, 'test', background)
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

    // profile --catalog on the same test reads, with two catalogs: A, every sequence of the release is
    // catalogued (upper bound); B, only what the benchmark does not score (no test read's own gene, but
    // not a floor: conserved regions still settle reads). Both are indexed as the release's catalog
    // is (ext.args); A is also run as profile --catalog in full on a few samples.
    BENCHMARK_CATALOGGENES(refs, families, BENCHMARK_SIMULATE_TEST.out.reads)
    catalog_sets = BENCHMARK_CATALOGGENES.out.genes_a.mix(BENCHMARK_CATALOGGENES.out.genes_b)
        .flatMap { scen, genes ->
            scen == 'A'
                ? [tuple([id: 'A', verify: ['opuAC.negative', 'proX.positive', 'trkH.negative']], genes)]
                : [tuple([id: 'B', verify: []], genes)]
        }
    OSMOTOOL_BUILDCATALOG_BENCHMARK(catalog_sets, db_calib)
    BENCHMARK_CATALOGMAP(BENCHMARK_SIMULATE_TEST.out.reads, db, OSMOTOOL_BUILDCATALOG_BENCHMARK.out.catalog)
    BENCHMARK_CATALOGSCORE(BENCHMARK_SIMULATE_TEST.out.reads, BENCHMARK_READS_TEST.out.results, cascade, families,
        BENCHMARK_CATALOGMAP.out.runs.map { _meta, d -> d }.collect())

    emit:
    db            = db                            // the release as osmotool reads it: [files..., hmms/]
    db_annotate   = db_calib                      // the same without the cascade: what annotation and build-catalog need,
                                                  // ready before the read benchmark
    fetch_manifest = REFS_BUILD.out.manifest      // refs/manifest.tsv: the UniProt release and fetch dates
    refs          = refs
    alignments    = HMM_BUILD.out.alignments
    hmms          = hmms
    release_files = release_files                 // the DIAMOND database and its files
    cascade       = cascade                       // cascade config and report levels
    done          = BENCHMARK_METRICS_TEST.out.metrics.mix(BENCHMARK_ANNOTATESCORE.out.annotate, BENCHMARK_CATALOGSCORE.out.metrics).collect()
}
