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

    // a release is published once: a second run of the same name only with -resume (after the gate)
    published = file("${params.repo}/releases/${params.release}")
    if( published.exists() && !workflow.resume )
        error "${published} exists: name a new release (--release), or add -resume to continue that run"

    emit:
    families      = channel.value(file(params.families_file, checkIfExists: true))
    source        = channel.value(source)
    pool          = channel.value(random_panel.resolve('pool.hmm'))
    pool_families = channel.value(random_panel.resolve('pool_families.tsv'))
    snapshot      = channel.value(snapshot)
}
