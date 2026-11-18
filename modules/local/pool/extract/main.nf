// The random-panel pool's models: the families of assets/random_panel/pool_families.tsv (the pool the
// BacDive analysis was specified with: 174 non-osmotic Pfam families, drawn from Pfam-A 39.0 and
// screened on the 43 pilot proteomes), taken from Pfam-A as they are, in Pfam's order.
// pool.hmm is what BACDIVE_PROTEINS searches (BACDIVE_ANALYZE draws its panels from pool_families.tsv);
// published to <repo>/<data_dir>/random_panel (conf/modules.config).
process POOL_EXTRACT {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    path(pfam)                                              // POOL_PFAM: Pfam-A.hmm.gz
    path(families, stageAs: 'pool_families.tsv')            // the pool's families

    output:
    path('pool.hmm'),          emit: hmm
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python

    script:
    """
    python3 <<'PYEOF'
    import csv
    import gzip

    with open("pool_families.tsv") as fh:
        wanted = {r["pfam"] for r in csv.DictReader(fh, delimiter="\\t")}
    found = set()
    with gzip.open("${pfam}", "rt") as fin, open("pool.hmm", "w") as out:
        record, acc = [], None
        for line in fin:
            record.append(line)
            if line.startswith("ACC "):
                acc = line.split(None, 1)[1].strip().split(".")[0]
            elif line.startswith("//"):
                if acc in wanted:
                    out.writelines(record)
                    found.add(acc)
                record, acc = [], None
    missing = sorted(wanted - found)
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} pool families not in ${pfam}: {missing[:10]}")
    print(f"pool: {len(found)} families -> pool.hmm")
    PYEOF
    """

    stub:
    """
    printf 'HMMER3/f\\nNAME  stub\\n//\\n' > pool.hmm
    """
}
