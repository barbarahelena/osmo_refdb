// The catalog benchmark's catalog decisions: every test sample mapped to one benchmark catalog with
// profile --catalog --catalog_only (no translated search, no cascade), so each settled read's
// decision is in <sample>.catalog_hits.tsv. With meta.verify, also the full profile --catalog
// (translated fallback, the release's cascade) on a few samples, which BENCHMARK_CATALOGSCORE checks
// its per-read rebuild against. ext.args: the read benchmark's DIAMOND filters, for the verify runs.
// The database decides whether catalog calls fall through: BENCHMARK_CATALOGMAP gets the release's db,
// with profile_report_levels.tsv, so a unit with catalog_reported no has its reads left out of
// catalog_hits.tsv (osmotool >= 0.11.0). BENCHMARK_CATALOGMAP_CALIB gets
// db_calib, made before the levels exist, so it sees every catalog call the decision is made on.
process BENCHMARK_CATALOGMAP {
    // calls: osmotool profile (bwa-mem3, diamond, hmmsearch)
    container 'docker.io/barbarahelena/osmotool:0.11.1'
    label 'process_medium'
    tag "${meta.id}"
    input:
    path(reads)
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE), staged as db/
    tuple val(meta), path(catalog, stageAs: 'catalog')

    output:
    tuple val(meta), path("${meta.id}"), emit: runs   // <id>/runs/<sample>.*, and <id>/verify/<sample>.* with meta.verify
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool

    script:
    def args   = task.ext.args ?: ''
    def verify = (meta.verify ?: []).join(' ')
    """
    mkdir -p ${meta.id}/runs ${meta.id}/verify
    for r1 in ${reads}/*_R1.fastq; do
        s=\$(basename "\$r1" _R1.fastq)
        osmotool profile db -1 "\$r1" -2 ${reads}/\${s}_R2.fastq --out_prefix ${meta.id}/runs/\$s \\
            --catalog_dir catalog --catalog_only --no_cascade --no_markers --keep_aln --threads ${task.cpus}
    done
    # the verify samples named in meta.verify that exist; with a families subset none may, then the first three
    todo=""
    for s in ${verify}; do [ -f ${reads}/\${s}_R1.fastq ] && todo="\$todo \$s"; done
    if [ -n "${verify}" ] && [ -z "\$todo" ]; then
        todo=\$(ls ${reads}/*_R1.fastq | head -3 | xargs -n1 basename | sed 's/_R1.fastq\$//')
    fi
    for s in \$todo; do
        osmotool profile db -1 ${reads}/\${s}_R1.fastq -2 ${reads}/\${s}_R2.fastq --out_prefix ${meta.id}/verify/\$s \\
            --catalog_dir catalog --min_identity 0.40 --min_query_cover 0.50 --no_markers --threads ${task.cpus} ${args}
    done
    """

    stub:
    """
    mkdir -p ${meta.id}/runs ${meta.id}/verify
    """
}
