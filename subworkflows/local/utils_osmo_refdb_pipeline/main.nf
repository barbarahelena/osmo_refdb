//
// Parameter checks before anything is submitted, and the helpers they use.
//

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

// The sets a family's model is scored on for its calibration, from its refs files by kind
// ([kind: file], kind = 'positive.calib', 'negative.train', ...; empty files left out), in the order
// the HMM_BUILD subworkflow used: positives from calib, else test (the benchmark's positives overlap), else all;
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

// A refs file's family and kind: <family>.<kind>.faa (ectA.positive.calib.faa -> ectA, positive.calib).
def refs_kind(f) {
    def base = f.name - ~/\.faa$/
    def i = base.indexOf('.')
    i < 0 ? tuple(base, '', f) : tuple(base.substring(0, i), base.substring(i + 1), f)
}

// The sets a family's DIAMOND cutoff is calibrated on, as calibrate_diamond_cutoffs.py picks them
// from the family's refs files by kind ([kind: file], empty files included): positive.calib if it exists,
// else positive.test; negative.train if it exists, else all negatives. An empty or missing set is not
// searched (08b scores it as no sequences).
def diamond_sets(sets) {
    def pos = sets['positive.calib'] ?: sets['positive.test']
    def neg = sets['negative.train'] ?: sets['negative']
    [positive: pos, negative: neg].findAll { _kind, f -> f && f.size() > 0 }
}

// One family's alignment: 'skip' with fewer than 2 train sequences (no alignment, and so no hmmbuild),
// else 'align'.
def alignment_plan(train) {
    train.readLines().count { line -> line.startsWith('>') } < 2 ? 'skip' : 'align'
}

def check_catalog_name() {
    if( file(params.catalog_genes).name != 'GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz' )
        error "--catalog_genes must end in GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz (that is what the download is called)"
}

workflow PIPELINE_INITIALISATION {
    main:
    require_params(['repo', 'release', 'catalog_name', 'catalog_genes', 'full_index_dir', 'data_dir', 'nchunks'])
    check_catalog_name()
    // harvested BacDive batches to start from (--bacdive_full), or []
    snapshot = params.bacdive_snapshot ? file(params.bacdive_snapshot, checkIfExists: true) : []

    emit:
    families = channel.value(file(params.families_file, checkIfExists: true))
    markers  = channel.value(file(params.markers_file, checkIfExists: true))
    snapshot = channel.value(snapshot)
}
