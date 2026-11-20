//
// Real-genome validation of a release: does it call genes correctly on real, independently annotated
// genomes? Two checks, both on inputs committed in the repo, so it runs for every release:
//   1. system calls on the reference genomes (--reference_genomes_dir, osmotool's test_genomes/)
//      against assets/real_genome_validation/expected_systems.tsv
//   2. read-level recall and precision on reads simulated from the genomes with curated target genes
//      (assets/real_genome_validation/target_orfs.tsv, genomes in --test_genomes_dir), with translated
//      search only and with the release's catalog (the index it ships: what profile uses by default)
//
include { OSMOTOOL_ANNOTATE as OSMOTOOL_ANNOTATE_GENOME } from '../../../modules/local/osmotool/annotate/main'
include { OSMOTOOL_PROFILE  as OSMOTOOL_PROFILE_RECALL  } from '../../../modules/local/osmotool/profile/main'
include { OSMOTOOL_PROFILE  as OSMOTOOL_PROFILE_RECALL_CATALOG } from '../../../modules/local/osmotool/profile/main'
include { WGSIM                   } from '../../../modules/local/wgsim/main'
include { ORFM                    } from '../../../modules/local/orfm/main'
include { HMMER_HMMSCAN           } from '../../../modules/local/hmmer/hmmscan/main'
include { VALIDATION_CHECKSYSTEMS } from '../../../modules/local/validation/checksystems/main'
include { VALIDATION_TARGETCOORDS } from '../../../modules/local/validation/targetcoords/main'
include { VALIDATION_READRECALL   } from '../../../modules/local/validation/readrecall/main'

workflow REAL_GENOME_VALIDATION {
    take:
    db          // the release as osmotool reads it
    hmms        // HMM_CALIBRATE: the pressed HMM database, for hmmscan on the reads
    catalog     // CATALOG_PACKAGE: [meta, catalogs/<name>/]

    main:
    targets  = file("${projectDir}/assets/real_genome_validation/target_orfs.tsv", checkIfExists: true)
    expected = file("${projectDir}/assets/real_genome_validation/expected_systems.tsv", checkIfExists: true)
    test_dir = params.test_genomes_dir ?: "${projectDir}/assets/test_genomes"

    reference = params.reference_genomes_dir
        ? channel.fromPath("${params.reference_genomes_dir}/*.fna").map { f -> tuple([id: f.baseName], f) }
        : channel.empty()
    curated = channel.fromPath(targets)
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.genome }
        .unique()
        .map { g -> tuple([id: g], file("${test_dir}/${g}.fna", checkIfExists: true)) }

    // one annotation per genome, whichever check needs it (Prodigal proteins kept, ext.args)
    genomes = reference.mix(curated).unique { pair -> pair[0].id }
    OSMOTOOL_ANNOTATE_GENOME(genomes.map { meta, fna -> tuple(meta, fna, []) }, db)

    // 1. system calls
    VALIDATION_CHECKSYSTEMS(
        OSMOTOOL_ANNOTATE_GENOME.out.systems.mix(OSMOTOOL_ANNOTATE_GENOME.out.gene_counts)
            .map { _meta, f -> f }.collect(),
        expected)

    // 2. read recall, on the genomes that have curated target genes
    curated_proteins = OSMOTOOL_ANNOTATE_GENOME.out.proteins.join(curated)   // meta, faa, fna
    VALIDATION_TARGETCOORDS(curated_proteins.map { meta, faa, _fna -> tuple(meta, faa) }, targets)
    WGSIM(curated_proteins.map { meta, _faa, fna -> tuple(meta, fna) })
    reads = WGSIM.out.reads

    OSMOTOOL_PROFILE_RECALL(reads, db)   // DIAMOND, alignments kept (ext.args)
    // the same with the catalog, staged in the database as db/<name>/ (ext.args: --catalog db/<name>)
    OSMOTOOL_PROFILE_RECALL_CATALOG(reads, db.mix(catalog.map { _meta, d -> d }).collect())
    catalog_calls = OSMOTOOL_PROFILE_RECALL_CATALOG.out.results
        .map { meta, files -> tuple(meta, files.find { f -> f.name.endsWith('.catalog_hits.tsv') } ?: [],
                                          files.find { f -> f.name.endsWith('.blastx.tsv') } ?: []) }
    ORFM(reads)
    HMMER_HMMSCAN(ORFM.out.orfs, hmms, 'osmo_refdb.hmm')   // HMM on the same reads

    VALIDATION_READRECALL(
        reads.join(VALIDATION_TARGETCOORDS.out.regions)
             .join(OSMOTOOL_PROFILE_RECALL.out.alignments)
             .join(HMMER_HMMSCAN.out.tblout)
             .join(catalog_calls))

    emit:
    report   = VALIDATION_CHECKSYSTEMS.out.report
                   .mix(VALIDATION_READRECALL.out.results.map { _meta, truth, cmp -> [truth, cmp] }.flatten())
                   .collect()
}
