//
// Parameter checks before anything is submitted, and the helpers they use.
//

def under_repo(p) {
    p.toString().startsWith('/') ? p.toString() : "${params.repo}/${p}".toString()
}

def require_params(names) {
    names.each { n ->
        if( !params[n] )
            error "--${n} is required (set it in nextflow.config, or on the command line for --release)"
    }
}

// The families in families.yaml, as one meta map each: [id, trim_gt (trimAl's -gt, default 0.8),
// pfam_model (a curated Pfam model used in place of hmmbuild, or '')].
def family_metas(families_yaml) {
    def data = new org.yaml.snakeyaml.Yaml().load(families_yaml.text)
    data.families.collect { fam -> [id: fam.name, trim_gt: fam.trim_gt ?: 0.8, pfam_model: fam.pfam_model ?: ''] }
}

// The HMM decoys in families.yaml (hmm_decoys): [id, trim_gt, competes_with, decoy_aware (decoy_aware_ga)].
def decoy_metas(families_yaml) {
    def data = new org.yaml.snakeyaml.Yaml().load(families_yaml.text)
    (data.hmm_decoys ?: []).collect { d ->
        [id: d.name, trim_gt: d.trim_gt ?: 0.8, competes_with: d.competes_with ?: [], decoy_aware: d.decoy_aware_ga ?: false]
    }
}

def sha256(path) {
    java.security.MessageDigest.getInstance('SHA-256').digest(path.bytes).encodeHex().toString()
}

// The sets a family's model is scored on for its calibration, from its refs files by kind
// ([kind: file], kind = 'positive.calib', 'negative.train', ...; empty files left out), in the order
// 05_build_hmms.sh used: positives from calib, else test (the benchmark's positives overlap), else all;
// negatives from hmm_negative (an HMM-only pool), else negative.train, else all negatives.
def scoring_sets(family, sets) {
    def pos = ['positive.calib', 'positive.test', 'positive'].find { k -> sets.containsKey(k) }
    def neg = ['hmm_negative', 'negative.train', 'negative'].find { k -> sets.containsKey(k) }
    if( pos && pos != 'positive.calib' )
        log.warn "[${family}] no positive.calib split: calibrating on ${pos} (benchmark positives overlap)"
    if( neg == 'negative' )
        log.warn "[${family}] no negative.train split: calibrating on the full negative set (benchmark negatives overlap)"
    [positive: pos ? sets[pos] : null, negative: neg ? sets[neg] : null]
}

def source_manifest_row(source, id) {
    def manifest = source.resolve('alignments/manifest.tsv')
    if( !manifest.exists() )
        return null
    manifest.readLines().drop(1).collect { line -> line.split('\t') }.find { cols -> cols[0] == id }
}

// What to do with one family's alignment: 'skip' (fewer than 2 sequences: no alignment, and so no
// hmmbuild), 'reuse' the source release's (multi-threaded MAFFT is not reproducible run to run, so an
// unchanged family keeps the alignment it was built with), or 'align'. Reused only when the train set
// is byte-identical to the source's, both of the source's alignment files exist, and the source's
// alignments/manifest.tsv, where it has one, records the same trim_gt.
def alignment_plan(meta, train, source) {
    if( train.readLines().count { line -> line.startsWith('>') } < 2 )
        return 'skip'
    def old_train = source.resolve("refs/${meta.id}.positive.train.faa")
    def old_files = ['aln.fasta', 'trimmed.fasta'].collect { ext -> source.resolve("alignments/${meta.id}.${ext}") }
    if( !old_train.exists() || !old_files.every { f -> f.exists() } || sha256(old_train) != sha256(train) )
        return 'align'
    def row = source_manifest_row(source, meta.id)
    if( source.resolve('alignments/manifest.tsv').exists() && (!row || row[2] != meta.trim_gt.toString()) )
        return 'align'
    return 'reuse'
}

// The source release's files for an HMM decoy, if it has all of them: its input set, its clustering
// and its alignment (named <decoy>.aln / .trimmed.aln before release v11, .aln.fasta / .trimmed.fasta since).
def first_existing(source, names) {
    names.collect { n -> source.resolve(n) }.find { f -> f.exists() }
}

def source_decoy_files(source, id) {
    [input:     source.resolve("refs/${id}.hmm_decoy.faa"),
     clustered: source.resolve("refs/${id}.hmm_decoy.clustered.faa"),
     clusters:  source.resolve("refs/${id}.hmm_decoy.clustered.faa.clstr"),
     raw:       first_existing(source, ["alignments/${id}.aln.fasta", "alignments/${id}.aln"]),
     trimmed:   first_existing(source, ["alignments/${id}.trimmed.fasta", "alignments/${id}.trimmed.aln"])]
}

// An HMM decoy is rebuilt from scratch ('align': clustering, alignment, trimming) unless its input
// set is byte-identical to the source's and the source has its clustering and alignment ('reuse'),
// with the same trim_gt where the source's manifest records it.
def decoy_plan(meta, input, source) {
    def old = source_decoy_files(source, meta.id)
    if( !old.input.exists() || !old.clustered.exists() || !old.raw || !old.trimmed || sha256(old.input) != sha256(input) )
        return 'align'
    def row = source_manifest_row(source, meta.id)
    if( row && row[2] != meta.trim_gt.toString() )
        return 'align'
    return 'reuse'
}

def check_catalog_name() {
    if( file(params.catalog_genes).name != 'GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz' )
        error "--catalog_genes must end in GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz (that is what the download is called)"
}

workflow PIPELINE_INITIALISATION {
    main:
    require_params(['repo', 'release', 'reuse_refs_from', 'catalog_name',
                    'catalog_genes', 'full_index_dir', 'data_dir', 'nchunks'])
    check_catalog_name()

    // the inputs that are not made by the pipeline: snapshots, made once by hand or shipped
    source = file(under_repo(params.reuse_refs_from))
    if( !source.resolve('refs/manifest.tsv').exists() )
        error "no release snapshot at ${source} (refs/manifest.tsv missing; rsync it, see docs/hpc.md)"
    random_panel = file(params.random_panel ?: under_repo("${params.data_dir}/random_panel"))
    ['pool.hmm', 'pool_families.tsv'].each { f ->
        if( !random_panel.resolve(f).exists() )
            error "no random-panel pool snapshot: ${random_panel}/${f} missing (--random_panel, see docs/hpc.md)"
    }
    snapshot = params.bacdive_snapshot ? file(params.bacdive_snapshot, checkIfExists: true) : []

    emit:
    families      = channel.value(file(params.families_file, checkIfExists: true))
    source        = channel.value(source)
    pool          = channel.value(random_panel.resolve('pool.hmm'))
    pool_families = channel.value(random_panel.resolve('pool_families.tsv'))
    snapshot      = channel.value(snapshot)
}
