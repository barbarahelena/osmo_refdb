//
// A release's references, from scratch (UniProt and RefSeq are live: what they hold that day; the
// release and fetch dates are in refs/manifest.tsv):
//   REFS_FETCHFAMILY (per family) + REFS_FETCHSHARED (fusion pairs, HMM decoys) -> REFS_COMBINE -> REFS_PURITY -> REFS_LENGTH -> REFS_REFSEQ -> REFS_BACKGROUND (when asked for)
//   -> REFS_CLUSTER, per family's positives and negatives
// Each step works on a copy of the previous one's refs/. For families.yaml and for markers.yaml.
//
include { UNIPROT_RELEASE } from '../../../modules/local/uniprot/release/main'
include { REFS_FETCHFAMILY } from '../../../modules/local/refs/fetchfamily/main'
include { REFS_FETCHSHARED } from '../../../modules/local/refs/fetchshared/main'
include { REFS_COMBINE     } from '../../../modules/local/refs/combine/main'
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
    // Each family fetched on its own, keyed for --refs_store on the UniProt release, that family's
    // entry and the fetch script: a changed family refetches only itself. The fusion-pair candidates and
    // the HMM decoys need every family's fetch, keyed on the decoys, the fusion_partner families and the script.
    UNIPROT_RELEASE()
    release = UNIPROT_RELEASE.out.release.map { r -> r.trim() }
    script_md5 = file("${projectDir}/modules/local/refs/templates/fetch_refs.py").text.md5()
    entries = families.flatMap { f -> family_entries(f) }
    REFS_FETCHFAMILY(release.combine(entries)
        .map { r, name, entry -> tuple(r, name, (entry + script_md5).md5().take(16), entry) })
    fetched = REFS_FETCHFAMILY.out.refs.collect(sort: true)
    REFS_FETCHSHARED(release.combine(families).map { r, f -> tuple(r, shared_key(f, script_md5), f) }, fetched)
    REFS_COMBINE(release, families, fetched, REFS_FETCHSHARED.out.refs)
    REFS_PURITY(REFS_COMBINE.out.refs, families)
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

// Each family of a families file as [name, its own one-family families file (YAML)].
def family_entries(families_yaml) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def data = yaml.load(families_yaml.text)
    data.families.collect { fam -> [fam.name, yaml.dump([families: [fam]])] }
}

// The shared fetch's store key: the HMM decoys, the fusion_partner families' entries, the fetch script.
def shared_key(families_yaml, script_md5) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def data = yaml.load(families_yaml.text)
    def fused = data.families.findAll { fam -> fam.fusion_partner }
    (yaml.dump([hmm_decoys: data.hmm_decoys ?: [], fused: fused]) + script_md5).md5().take(16)
}
