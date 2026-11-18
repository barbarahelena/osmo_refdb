//
// A release's references, from scratch (UniProt and RefSeq are live: what they hold that day; the
// release and fetch dates are in refs/manifest.tsv):
//   REFS_FETCH -> REFS_PURITY -> REFS_LENGTH -> REFS_REFSEQ -> REFS_BACKGROUND (when asked for)
//   -> REFS_CLUSTER, per family's positives and negatives
// Each step works on a copy of the previous one's refs/. For families.yaml and for markers.yaml.
//
include { UNIPROT_RELEASE } from '../../../modules/local/uniprot/release/main'
include { REFS_FETCH      } from '../../../modules/local/refs/fetch/main'
include { REFS_PURITY     } from '../../../modules/local/refs/purity/main'
include { REFS_LENGTH     } from '../../../modules/local/refs/length/main'
include { REFS_REFSEQ     } from '../../../modules/local/refs/refseq/main'
include { REFS_BACKGROUND } from '../../../modules/local/refs/background/main'
include { REFS_CLUSTER    } from '../../../modules/local/refs/cluster/main'

workflow REFS_BUILD {
    take:
    families        // families.yaml or markers.yaml
    background      // true: also the benchmark's random background (the markers have no benchmark)

    main:
    // keyed for --refs_store on the UniProt release, the families file and the fetch script
    UNIPROT_RELEASE()
    fetch_script = file("${projectDir}/modules/local/refs/templates/fetch_refs.py")
    REFS_FETCH(UNIPROT_RELEASE.out.release.map { r -> r.trim() }.combine(families)
        .map { r, f -> tuple(r, (f.text + fetch_script.text).md5().take(16), f) })
    REFS_PURITY(REFS_FETCH.out.refs, families)
    REFS_LENGTH(REFS_PURITY.out.refs, families)
    REFS_REFSEQ(REFS_LENGTH.out.refs, families)
    if( background ) {
        REFS_BACKGROUND(REFS_REFSEQ.out.refs, families)
        files = REFS_BACKGROUND.out.files.flatten()
    } else {
        files = REFS_REFSEQ.out.files.flatten()
    }

    // each family's positives and negatives clustered; the rest as they are
    sets = files.filter { f -> f.name ==~ /[^.]+\.(positive|negative)\.faa/ && f.size() > 0 }
    REFS_CLUSTER(sets.map { f -> tuple([id: f.name - '.faa'], f) })
    refs = files.filter { f -> !(f.name ==~ /[^.]+\.(positive|negative)\.faa/ && f.size() > 0) }
        .mix(REFS_CLUSTER.out.clustered.flatMap { _meta, faa, clustered -> [faa, clustered] },
             REFS_CLUSTER.out.clusters.map { _meta, clstr -> clstr })
        .collect()

    emit:
    refs     = refs                                                       // every refs/ file
    manifest = files.filter { f -> f.name == 'manifest.tsv' }.first()    // the fetch's: UniProt release, dates
}
