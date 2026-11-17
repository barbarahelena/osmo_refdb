#!/usr/bin/env nextflow
/*
 * osmo_refdb on a SLURM cluster, as one workflow. Every step is a job in a container of its own
 * (declared in its module; containers/README.md says what each needs). Start it from a SLURM job of
 * your own (a small one, it only submits and watches the others), not on a login node.
 *
 *   downloads: GMGC human gut | BacDive strains + genomes
 *   then, in parallel:
 *     full catalog index                       (needs only the gene file)
 *     BacDive proteins x nchunks               (needs only the genomes)
 *     release build -> catalog annotation -> screen            <- THE GATE
 *   after the gate (--gate_passed true):
 *     BacDive annotation x nchunks -> analysis + NaCl model
 *   last: PROVENANCE.json
 *
 * The gate is a person reading the screen, so the workflow stops there on the first run. Run the
 * same command again with `-resume --gate_passed true`: the cached steps are skipped and the
 * BacDive part runs on the frozen release.
 *
 * Other entry:  -entry profile_samples --samples s.csv    osmotool profile, one job per sample
 *               (CSV columns id,r1,r2; results in <outdir>/profile)
 *
 * Layout follows the nf-core template: workflows/, subworkflows/local/, modules/local/<step>/,
 * conf/base.config (resources), conf/modules.config (per-step settings). What it does not use: the
 * nf-schema plugin (nextflow_schema.json documents the parameters but is not enforced), conda
 * (osmotool is only distributed as a container), nf-core utility subworkflows, nf-test.
 *
 * Skipping finished work: -resume reuses results from the work directory; storeDir keeps a task's
 * output in a permanent place and skips the task whenever it is already there, even without work/.
 * storeDir is used where the output is a real file: the gene file and the full catalog index.
 * Everything else (BacDive data, the release, the catalog annotation) is modified in place by shell
 * scripts, so those tasks only pass a name on and check their product exists, and rely on -resume.
 * Their scripts skip work already done, so rerunning after deleting work/ is cheap for BacDive. If you
 * delete a release, delete the work directory (or run without -resume) too. See docs/hpc.md.
 *
 * Nextflow starts Apptainer with --no-home and `env -`: neither the host home directory nor the host
 * environment is inside the container.
 */

include { PIPELINE_INITIALISATION; require_params } from './subworkflows/local/utils_osmo_refdb_pipeline/main'
include { OSMO_REFDB      } from './workflows/osmo_refdb'
include { PROFILE_SAMPLE  } from './modules/local/profile_sample/main'

workflow {
    main:
    PIPELINE_INITIALISATION()
    OSMO_REFDB()
}

workflow profile_samples {
    main:
    require_params(['repo', 'release', 'catalog_name', 'full_index_dir', 'samples'])
    samples = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .map { row -> tuple([id: row.id], file(row.r1, checkIfExists: true), file(row.r2, checkIfExists: true)) }
    PROFILE_SAMPLE(samples)
}
