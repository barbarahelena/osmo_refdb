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

// The families in families.yaml, as one meta map each: [id, trim_gt (trimAl's -gt, default 0.8)].
def family_metas(families_yaml) {
    def data = new org.yaml.snakeyaml.Yaml().load(families_yaml.text)
    data.families.collect { fam -> [id: fam.name, trim_gt: fam.trim_gt ?: 0.8] }
}

def sha256(path) {
    java.security.MessageDigest.getInstance('SHA-256').digest(path.bytes).encodeHex().toString()
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
    def manifest = source.resolve('alignments/manifest.tsv')
    if( manifest.exists() ) {
        def row = manifest.readLines().drop(1).collect { line -> line.split('\t') }.find { cols -> cols[0] == meta.id }
        if( !row || row[2] != meta.trim_gt.toString() )
            return 'align'
    }
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
