// BacDive strains and genomes: harvest the BacDive API, the phenotype table, the genomes from NCBI,
// and the within-domain salt-in control (its strains, and genomes for those the main set lacks).
// snapshot: a bacdive_raw/ directory of harvested batches (--bacdive_snapshot); batches in it are not
// fetched again, so a shipped snapshot reproduces the published table (BacDive itself is live). [] to
// harvest everything. The tables and the batches are published to <repo>/<data_dir> (conf/modules.config).
process BACDIVE_FETCH {
    label 'process_high'
    errorStrategy 'retry'   // a dropped connection is not a failure
    maxRetries 3

    // calls: curl, python (stdlib only), awk, gzip, xargs, sed, grep
    // Seqera container 'bacdive-fetch' (environment.yml next to this file: its packages)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/curl_gawk_python:2d801c4a5165bf61' :
        'community.wave.seqera.io/library/curl_gawk_python:de4685943c4289cd' }"

    input:
    path(snapshot, stageAs: 'snapshot')

    output:
    path('bacdive/genomes'),                              emit: genomes
    path('bacdive/pilot_strains.tsv'),                    emit: pilot
    path('bacdive/bacdive_classes.tsv'),                  emit: classes
    path('bacdive/strategy_control/genomes'),             emit: control_genomes
    path('bacdive/strategy_control/pilot_strains.tsv'),   emit: control_pilot
    path('bacdive/bacdive_raw'),                          emit: batches
    path('bacdive/genome_fetch.log'), optional: true,     emit: log
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('curl'), eval("curl --version | head -1 | cut -d' ' -f2"), topic: versions, emit: versions_curl

    script:
    def S = "${moduleDir}/templates/phenotype_validation"
    """
    # Each strain's assembly from NCBI into <dir>/genomes/<bid>_<accession>.fna.gz. BacDive's INSDC
    # accessions are unversioned (GCA_013177455): the versioned FTP folder is read from the
    # accession's parent folder, newest version first, falling back to the next where a suppressed
    # version has no FASTA. Each download is checked (gzip -t) under a temporary name, so an interrupted
    # one never looks complete; genomes already there are kept. Logged to <dir>/genome_fetch.log.
    fetch_genomes() (   # <dir with pilot_strains.tsv> [bacdive_id...]; PAR downloads at once
    set +e -uo pipefail
    PV=\${1:?usage: fetch_genomes <outdir> [bacdive_id...]}; shift || true
    PAR=\${PAR:-4}
    TSV="\$PV/pilot_strains.tsv"
    [ -s "\$TSV" ] || { echo "missing \$TSV -- run build_table.py first" >&2; exit 1; }
    mkdir -p "\$PV/genomes"

    # look columns up by name: positions shift whenever build_table.py gains one
    # (the pilot's hard-coded \$9 silently became nogrowth_at_pct)
    pairs=\$(awk -F'\\t' '
      NR==1 { for (i=1;i<=NF;i++) col[\$i]=i
              if (!col["bacdive_id"] || !col["insdc"]) { print "ERR" > "/dev/stderr"; exit 1 }
              next }
      { print \$col["bacdive_id"] "\\t" \$col["insdc"] }' "\$TSV") \\
      || { echo "\$TSV lacks a bacdive_id or insdc column" >&2; exit 1; }
    if [ "\$#" -gt 0 ]; then
      want=" \$* "
      pairs=\$(printf '%s\\n' "\$pairs" | awk -F'\\t' -v w="\$want" 'index(w, " " \$1 " ")')
    fi

    one() {
      bid="\$1"; acc="\$2"
      out="\$G/\${bid}_\${acc}.fna.gz"
      { [ -s "\$out" ] || [ -s "\$G/\${bid}_\${acc}.fna" ]; } && { echo "CACHED \$bid"; return 0; }
      p=\$(echo "\${acc#GCA_}" | sed 's/\\(...\\)\\(...\\)\\(...\\)/\\1\\/\\2\\/\\3/')
      base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/\$p/"
      dirs=\$(curl -sS -m 60 --retry 3 --retry-delay 2 "\$base" \\
             | grep -oE "\${acc}\\.[0-9]+_[^/\\"]+" | sort -u -V -r)
      [ -z "\$dirs" ] && { echo "FAIL \$bid \$acc: no FTP dir (suppressed?)"; return 0; }
      tmp="\$out.part"
      # newest version first; a suppressed version keeps its folder but only
      # the metadata reports (the FASTA 404s), so fall back to the next one
      for dir in \$dirs; do
        if curl -sS -f -m 600 --retry 3 --retry-delay 2 -o "\$tmp" \\
             "\${base}\${dir}/\${dir}_genomic.fna.gz" 2>/dev/null && gzip -t "\$tmp" 2>/dev/null; then
          mv "\$tmp" "\$out"; echo "OK \$bid \$dir"; return 0
        fi
        rm -f "\$tmp"
      done
      echo "FAIL \$bid \$acc: no genomic FASTA in any version (\$(echo \$dirs | tr ' ' ','))"
    }
    export -f one
    export G="\$PV/genomes"

    # each line is "<bacdive_id>\\t<accession>"; -L 1 hands both fields to one()
    LOG="\$PV/genome_fetch.log"
    printf '%s\\n' "\$pairs" | grep -v '^\$' \\
      | xargs -P "\$PAR" -L 1 bash -c 'one "\$1" "\$2"' _ | tee -a "\$LOG" \\
      | awk '{n[\$1]++} NR%250==0 {printf "  %d done (OK=%d CACHED=%d FAIL=%d)\\n", NR, n["OK"], n["CACHED"], n["FAIL"]; fflush()}
             \$1=="FAIL" {print; fflush()}
             END {printf "downloaded=%d cached=%d failed=%d\\n", n["OK"], n["CACHED"], n["FAIL"]}'
    )

    mkdir -p bacdive/bacdive_raw
    if [ -d snapshot ]; then cp -RL snapshot/. bacdive/bacdive_raw/; fi
    # 1. harvest (batches already there are skipped), 2. phenotype table, 3. genomes from NCBI
    bash ${S}/harvest_bacdive.sh bacdive sweep
    python3 ${S}/build_table.py --dir bacdive
    PAR=${task.cpus} fetch_genomes bacdive
    # 4. the within-domain salt-in control: strains, and genomes (links to the main set's where it has them)
    python3 ${S}/build_strategy_control.py --dir bacdive
    PAR=4 fetch_genomes bacdive/strategy_control
    test -d bacdive/genomes
    """

    stub:
    """
    mkdir -p bacdive/genomes bacdive/strategy_control/genomes bacdive/bacdive_raw
    echo '>c' > bacdive/genomes/1_GCA_000001.fna
    printf 'bacdive_id\\n1\\n' > bacdive/pilot_strains.tsv
    printf 'bacdive_id\\ttax_class\\n' > bacdive/bacdive_classes.tsv
    printf 'bacdive_id\\n' > bacdive/strategy_control/pilot_strains.tsv
    """
}
