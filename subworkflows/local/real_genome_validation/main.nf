//
// Real-genome validation of a release: does it call genes correctly on real, independently annotated
// genomes? Two checks, both on inputs committed in the repo, so it runs for every release:
//   1. system calls on the reference genomes (--reference_genomes_dir, osmotool's test_genomes/)
//      against assets/real_genome_validation/expected_systems.tsv
//   2. read-level recall and precision on reads simulated from the genomes with curated target genes
//      (assets/real_genome_validation/target_orfs.tsv, genomes in --test_genomes_dir)
//
include { ANNOTATE_GENOME          } from '../../../modules/local/annotate_genome/main'
include { CHECK_REFERENCE_SYSTEMS  } from '../../../modules/local/check_reference_systems/main'
include { TARGET_COORDS            } from '../../../modules/local/target_coords/main'
include { SIMULATE_GENOME_READS    } from '../../../modules/local/simulate_genome_reads/main'
include { REAL_READ_RECALL         } from '../../../modules/local/real_read_recall/main'

workflow REAL_GENOME_VALIDATION {
    take:
    release     // the release's name, once it is final

    main:
    targets  = file("${projectDir}/assets/real_genome_validation/target_orfs.tsv", checkIfExists: true)
    test_dir = params.test_genomes_dir ?: "${params.repo}/test_genomes"

    reference = params.reference_genomes_dir
        ? channel.fromPath("${params.reference_genomes_dir}/*.fna").map { f -> tuple(f.baseName, f) }
        : channel.empty()
    curated = channel.fromPath(targets)
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.genome }
        .unique()
        .map { g -> tuple(g, file("${test_dir}/${g}.fna", checkIfExists: true)) }

    // one annotation per genome, whichever check needs it
    ANNOTATE_GENOME(release, reference.mix(curated).unique { pair -> pair[0] })
    annotated = ANNOTATE_GENOME.out.annotated

    // 1. system calls
    expected = file("${projectDir}/assets/real_genome_validation/expected_systems.tsv", checkIfExists: true)
    CHECK_REFERENCE_SYSTEMS(annotated.map { _name, files -> files }.flatten().collect(), expected)

    // 2. read recall, on the genomes that have curated target genes
    curated_annotated = annotated.join(curated.map { name, fna -> tuple(name, fna) })
    TARGET_COORDS(curated_annotated.map { name, files, _fna -> tuple(name, files) }, targets)
    SIMULATE_GENOME_READS(curated_annotated.map { name, _files, fna -> tuple(name, fna) })
    REAL_READ_RECALL(release, SIMULATE_GENOME_READS.out.reads.join(TARGET_COORDS.out.regions))

    emit:
    report = CHECK_REFERENCE_SYSTEMS.out.report.mix(REAL_READ_RECALL.out.results.map { _name, files -> files }).collect()
}
