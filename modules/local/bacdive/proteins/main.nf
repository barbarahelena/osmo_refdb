// Release-independent, per genome: size/GC, Prodigal proteins (gzipped, for BACDIVE_ANNOTATE) and the
// random-panel pool search, for the BacDive genomes and the salt-in control's. Task i of nchunks takes
// every nchunks-th genome; each strain's results are one directory, <bid>/.
// Prodigal and hmmsearch are osmotool's own: annotate --proteins reads these proteins in its place.
process BACDIVE_PROTEINS {
    tag "chunk ${chunk}"
    label 'process_medium'

    // calls: prodigal, hmmsearch, awk, find, xargs, sort, gzip
    container 'docker.io/barbarahelena/osmotool:0.9.0'

    input:
    val(chunk)                                          // i/K: every K-th genome from the i-th
    path(genomes)                                       // BACDIVE_GENOMES
    path(control_genomes, stageAs: 'control_genomes')   // BACDIVE_GENOMES: the salt-in control's, or []
    path(pool)                                          // the random-panel pool, pool.hmm

    output:
    tuple val(chunk), path('strains'), path('control_strains'), emit: chunk   // for BACDIVE_ANNOTATE's chunk
    path('strains/*'),         optional: true, emit: strains                    // per strain, for the analysis
    path('control_strains/*'), optional: true, emit: control_strains
    tuple val("${task.process}"), val('prodigal'), eval("prodigal -v 2>&1 | grep -oE 'V[0-9.]+' | tr -d V"), topic: versions, emit: versions_prodigal
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    """
    # Per genome, <par> at a time, one thread each, into <out_dir>/<bid>/:
    #   <bid>.genome_stats.tsv  genome_bp (A/C/G/T only), gc, n_contigs, total_proteins: the size/GC baselines
    #   <bid>.pool_counts.tsv   proteins hit per random-pool family (hmmsearch --cut_ga); written last,
    #                           so it marks the strain done
    #   <bid>.faa.gz            Prodigal's proteins (-p meta, as osmotool annotate runs it), for BACDIVE_ANNOTATE
    # One genome's failure is reported (rc=...) and does not stop the others. CHUNK=i/K: every K-th genome.
    pool_search() (   # <genomes_dir> <out_dir> <pool.hmm> <parallel>
    set +e -uo pipefail
    GENOMES=\${1:?usage: pool_search <genomes_dir> <out_dir> <pool.hmm> [par]}
    OUT=\${2:?missing out_dir}
    POOL=\${3:?missing pool.hmm}
    PAR=\${4:-8}
    mkdir -p "\$OUT"

    one() {
      g="\$1"; OUT="\$2"; POOL="\$3"
      b=\$(basename "\$g"); bid=\${b%%_*}
      d="\$OUT/\$bid"
      [ -s "\$d/\$bid.pool_counts.tsv" ] && { echo "\$bid cached"; return 0; }
      mkdir -p "\$d"
      t0=\$SECONDS
      tmp="\$d/\$bid.input.fna"
      case "\$g" in *.gz) gunzip -c "\$g" > "\$tmp" ;; *) cp "\$g" "\$tmp" ;; esac || {
        echo "\$bid rc=input_failed"; rm -f "\$tmp"; return 0; }

      # genome size counts A/C/G/T only (N runs and IUPAC codes excluded)
      stats=\$(awk 'BEGIN{n=0; gc=0; at=0}
           /^>/ {n++; next}
           { s=toupper(\$0); gc += gsub(/[GC]/, "", s); at += gsub(/[AT]/, "", s) }
           END{ printf "%d\\t%.5f\\t%d", gc+at, (gc+at)?gc/(gc+at):0, n }' "\$tmp")

      faa="\$d/\$bid.faa"
      if ! prodigal -i "\$tmp" -a "\$faa" -o /dev/null -p meta -q 2>"\$d/\$bid.log"; then
        echo "\$bid rc=prodigal_failed"; rm -f "\$tmp" "\$faa"; return 0
      fi
      rm -f "\$tmp"
      if ! hmmsearch --cut_ga --cpu 1 --noali --tblout "\$d/\$bid.tbl" -o /dev/null \\
           "\$POOL" "\$faa" 2>>"\$d/\$bid.log"; then
        echo "\$bid rc=hmmsearch_failed"; rm -f "\$faa"; return 0
      fi
      # proteins hit per family; a protein can hit several families, counted once each
      awk '!/^#/ {k=\$3 "\\t" \$1; if (!(k in seen)) {seen[k]=1; c[\$3]++}}
           END{print "family\\tn_proteins"; for (f in c) print f "\\t" c[f]}' \\
        "\$d/\$bid.tbl" > "\$d/\$bid.pool_counts.tmp"
      nprot=\$(grep -c '^>' "\$faa")
      printf 'genome_bp\\tgc\\tn_contigs\\ttotal_proteins\\n%s\\t%s\\n' "\$stats" "\$nprot" \\
        > "\$d/\$bid.genome_stats.tsv"
      gzip -f "\$faa"; rm -f "\$d/\$bid.tbl"
      # pool_counts.tsv is the "done" marker checked above, so it lands last
      mv "\$d/\$bid.pool_counts.tmp" "\$d/\$bid.pool_counts.tsv"
      echo "\$bid rc=0 wall=\$((SECONDS-t0))s proteins=\$nprot"
    }
    export -f one

    # CHUNK="i/K" (SLURM arrays): process only the genomes with (position - 1) % K == i,
    # so K array tasks split the sorted genome list between them without overlap.
    chunk_filter() {
      if [ -n "\${CHUNK:-}" ]; then awk -v i="\${CHUNK%/*}" -v k="\${CHUNK#*/}" '(NR - 1) % k == i'; else cat; fi
    }

    T0=\$SECONDS
    INPUTS="\$OUT/.inputs\${CHUNK:+.\${CHUNK%/*}}"
    find "\$GENOMES/" -maxdepth 1 \\( -name '*.fna' -o -name '*.fna.gz' \\) | sort | chunk_filter > "\$INPUTS"
    xargs -P "\$PAR" -I{} bash -c 'one "\$@"' _ {} "\$OUT" "\$POOL" < "\$INPUTS"
    N=\$(wc -l < "\$INPUTS" | tr -d ' '); W=\$((SECONDS-T0))
    echo "=== \${PAR}x1: \${W}s wall for \$N genomes ==="
    )

    export CHUNK=${chunk}   # i/K
    mkdir -p control_genomes   # none in pilot mode
    pool_search ${genomes} strains ${pool} ${task.cpus}
    pool_search control_genomes control_strains ${pool} ${task.cpus}
    """

    stub:
    """
    mkdir -p strains/${chunk.tokenize('/')[0]}00 control_strains
    touch strains/${chunk.tokenize('/')[0]}00/${chunk.tokenize('/')[0]}00.faa.gz
    """
}
