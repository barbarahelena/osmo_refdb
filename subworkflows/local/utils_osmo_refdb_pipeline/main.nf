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
    if( !file(under_repo("${params.reuse_refs_from}/refs/manifest.tsv")).exists() )
        error "no release snapshot at ${under_repo(params.reuse_refs_from)} (refs/manifest.tsv missing; rsync it, see docs/hpc.md)"
    if( !file(under_repo("${params.data_dir}/random_panel/pool.hmm")).exists() )
        error "no random-panel pool snapshot at ${under_repo(params.data_dir)}/random_panel/pool.hmm (see docs/hpc.md)"
}
