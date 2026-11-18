//
// The canonical run: everything that makes a release, in dependency order.
//
include { PROVENANCE         } from '../modules/local/provenance/main'
include { CATALOG_GENES      } from '../subworkflows/local/catalog_genes/main'
include { BACDIVE_DATA       } from '../subworkflows/local/bacdive_data/main'
include { RELEASE_CATALOG    } from '../subworkflows/local/release_catalog/main'
include { BACDIVE_ANALYSIS   } from '../subworkflows/local/bacdive_analysis/main'
include { REAL_GENOME_VALIDATION } from '../subworkflows/local/real_genome_validation/main'
include { MARKERS            } from '../subworkflows/local/markers/main'

workflow OSMO_REFDB {
    take:
    families        // families.yaml
    markers         // markers.yaml
    snapshot        // harvested BacDive batches to start from, or []

    main:
    // independent of the release: start at once, run alongside the build
    CATALOG_GENES()
    BACDIVE_DATA(snapshot)
    MARKERS(markers)

    // the release, then the catalog on it, then the gate
    RELEASE_CATALOG(families, CATALOG_GENES.out.genes)
    RELEASE_CATALOG.out.report.view { f -> f.text }
    // the release as osmotool reads it, now with its markers/
    db = RELEASE_CATALOG.out.db.mix(MARKERS.out.markers).collect()

    // after the gate: BacDive and the NaCl model on the frozen release, and the real-genome checks
    done = RELEASE_CATALOG.out.done.mix(MARKERS.out.markers).collect()
    if( params.gate_passed ) {
        data = BACDIVE_DATA.out
        BACDIVE_ANALYSIS(db, data.chunks, data.genomes, data.control_genomes, data.pilot, data.classes,
                         data.control_pilot, data.strains, data.control_strains, data.pool_families)
        REAL_GENOME_VALIDATION(db, RELEASE_CATALOG.out.hmms)
        done = done.mix(BACDIVE_ANALYSIS.out.done, REAL_GENOME_VALIDATION.out.report).collect()
    }

    // every step's tool versions, from the `versions` topic (modules: `eval(...), topic: versions`),
    // as one YAML file, {step: {tool: version}}; the chunked steps report identical ones, kept once
    software_versions = channel.topic('versions')
        .map { process, tool, version -> tuple(process.tokenize(':').last(), "${tool}: \"${version.trim()}\"".toString()) }
        .unique()
        .groupTuple(sort: true)
        .map { process, tools -> "\"${process}\":\n" + tools.collect { line -> "    ${line}" }.join('\n') }
        .collectFile(name: 'software_versions.yml', newLine: true, sort: true)

    // the receipt, once everything this invocation does to the release is done
    PROVENANCE(params.release, RELEASE_CATALOG.out.refs, RELEASE_CATALOG.out.alignments,
               RELEASE_CATALOG.out.hmms, RELEASE_CATALOG.out.release_files, RELEASE_CATALOG.out.cascade,
               RELEASE_CATALOG.out.catalog, software_versions, families, file(projectDir), done)

    emit:
    provenance = PROVENANCE.out.provenance
}
