#!/usr/bin/env nextflow
/*
 * osmo_refdb on a SLURM cluster, as one workflow. Every step is a job in a container of its own
 * (declared in its module, with its packages in the module's environment.yml). Start it from a SLURM job of
 * your own (a small one, it only submits and watches the others), not on a login node.
 *
 *   downloads: GMGC human gut | BacDive strains + genomes | UniProt and RefSeq references
 *   then, in parallel:
 *     BacDive proteins x nchunks               (needs only the genomes)
 *     markers                                  (markers.yaml)
 *     release build -> catalog annotation -> catalog screen
 *   on the release: BacDive annotation x nchunks -> analysis (+ NaCl model with --bacdive_full),
 *     the real-genome checks
 *   last: PROVENANCE.json
 *
 * The catalog screen is printed and kept in the release (catalogs/<name>/screen.tsv): read it before
 * publishing. A flagged family that is a real confuser class means a decoy and a rebuild.
 *
 * Other entry:  -entry profile_samples --samples s.csv    osmotool profile, one job per sample
 *               (CSV columns id,r1,r2; results in <outdir>/profile)
 *
 * Layout follows the nf-core template: workflows/, subworkflows/local/, modules/local/<tool>/<subtool>/
 * (osmotool/annotate, osmotool/profile, hmmer/hmmscan, ...; pipeline steps grouped as bacdive/, benchmark/,
 * validation/, catalog/),
 * conf/base.config (resources), conf/modules.config (per-step settings). What it does not use: the
 * nf-schema plugin (nextflow_schema.json documents the parameters but is not enforced), conda
 * (osmotool is only distributed as a container), nf-core utility subworkflows. Tests: tests/ (nf-test).
 *
 * Skipping finished work: -resume reuses results from the work directory; storeDir keeps a task's
 * output in a permanent place and skips the task whenever it is already there, even without work/.
 * storeDir is used for the gene file and, with --refs_store, the UniProt fetch. Everything else goes from step to step
 * as outputs and inputs; no step reads a published file. The release is published into
 * releases/<release> by the steps that make its parts (conf/modules.config). See README.md.
 *
 * Nextflow starts Apptainer with --no-home and `env -`: neither the host home directory nor the host
 * environment is inside the container.
 */

include { PIPELINE_INITIALISATION; require_params } from './subworkflows/local/utils_osmo_refdb_pipeline/main'
include { OSMO_REFDB      } from './workflows/osmo_refdb'
include { OSMOTOOL_PROFILE } from './modules/local/osmotool/profile/main'

workflow {
    main:
    PIPELINE_INITIALISATION()
    init = PIPELINE_INITIALISATION.out
    OSMO_REFDB(init.families, init.markers, init.snapshot)
}

workflow profile_samples {
    main:
    require_params(['repo', 'release', 'catalog_name', 'samples'])
    samples = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .map { row -> tuple([id: row.id], [file(row.r1, checkIfExists: true), file(row.r2, checkIfExists: true)]) }
    db    = files("${params.repo}/releases/${params.release}/*", checkIfExists: true)   // staged as osmotool's DATABASE
    // with the release's catalog and the index it ships (--catalog and --keep_aln: conf/modules.config)
    OSMOTOOL_PROFILE(samples, db)
}
