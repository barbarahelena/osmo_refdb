// osmotool annotate with the release's HMMs, per genome, from the proteins BACDIVE_PROTEINS cached (no
// Prodigal rerun), for the BacDive genomes and the salt-in control's. Task i takes the same genomes as
// BACDIVE_PROTEINS task i. Each strain's annotation is one directory, <bid>/.
process BACDIVE_ANNOTATE {
    tag "chunk ${chunk}"
    label 'process_low'

    // calls: osmotool annotate (hmmsearch, diamond)
    container 'docker.io/barbarahelena/osmotool:0.11.1'

    input:
    tuple val(chunk), path(strains), path(control_strains)   // BACDIVE_PROTEINS, the same chunk
    path(genomes)                                              // BACDIVE_GENOMES
    path(control_genomes, stageAs: 'control_genomes')
    path(db, stageAs: 'db/*')   // the release as osmotool reads it (DATABASE): its DIAMOND files, cascade and hmms/, staged as db/

    output:
    path('annotate/*'),         optional: true, emit: annotate
    path('control_annotate/*'), optional: true, emit: control_annotate
    tuple val("${task.process}"), val('osmotool'), eval("osmotool --version | cut -d' ' -f2"), topic: versions, emit: versions_osmotool
    tuple val("${task.process}"), val('diamond'), eval("diamond version | grep -oE '[0-9][0-9.]*' | head -1"), topic: versions, emit: versions_diamond
    tuple val("${task.process}"), val('hmmer'), eval("hmmsearch -h | grep -oE 'HMMER [0-9.]+' | cut -d' ' -f2"), topic: versions, emit: versions_hmmer

    script:
    """
    # osmotool annotate per genome, <par> at a time, one thread each, from the Prodigal proteins
    # BACDIVE_PROTEINS cached (<strains_dir>/<bid>/<bid>.faa.gz, Prodigal -p meta as annotate runs it),
    # passed with --proteins so Prodigal is not rerun (checked on 3 genomes against v8: gene_counts,
    # systems and gene_coordinates byte-identical to annotating the genome itself). One genome's failure
    # is reported (rc=...) and does not stop the others. CHUNK=i/K: every K-th genome from the i-th.
    annotate_genomes() (   # <genomes_dir> <strains_dir> <out_dir> <release> <parallel>
    set +e -uo pipefail
    GENOMES=\${1:?usage: annotate_genomes <genomes_dir> <strains_dir> <out_dir> <release> [par]}
    STRAINS=\${2:?missing strains_dir}
    OUT=\${3:?missing out_dir}
    REL=\${4:?missing release dir}
    PAR=\${5:-8}
    mkdir -p "\$OUT"

    one() {
      g="\$1"; STRAINS="\$2"; OUT="\$3"; REL="\$4"
      b=\$(basename "\$g"); bid=\${b%%_*}
      d="\$OUT/\$bid"
      [ -s "\$d/\$bid.gene_counts.tsv" ] && { echo "\$bid cached"; return 0; }
      prot="\$STRAINS/\$bid/\$bid.faa.gz"
      [ -s "\$prot" ] || { echo "\$bid rc=no_cached_proteins"; return 0; }
      mkdir -p "\$d"
      t0=\$SECONDS
      input="\$g"
      case "\$g" in
        *.gz) input="\$d/\$bid.input.fna"
              gunzip -c "\$g" > "\$input" || { echo "\$bid rc=gunzip_failed"; rm -f "\$input"; return 0; } ;;
      esac
      gunzip -c "\$prot" > "\$d/\$bid.input.faa" || { echo "\$bid rc=protein_gunzip_failed"; rm -f "\$d/\$bid.input.faa"; return 0; }
      osmotool annotate "\$REL" "\$input" --proteins "\$d/\$bid.input.faa" \\
        --out_prefix "\$d/\$bid" --threads 1 >"\$d/\$bid.log" 2>&1
      rc=\$?
      [ "\$input" != "\$g" ] && rm -f "\$input"
      rm -f "\$d/\$bid.input.faa"
      echo "\$bid rc=\$rc wall=\$((SECONDS-t0))s"
      [ \$rc -ne 0 ] && tail -5 "\$d/\$bid.log"
      return 0
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
    xargs -P "\$PAR" -I{} bash -c 'one "\$@"' _ {} "\$STRAINS" "\$OUT" "\$REL" < "\$INPUTS"
    N=\$(wc -l < "\$INPUTS" | tr -d ' '); W=\$((SECONDS-T0))
    echo "=== \${PAR}x1 (cached proteins): \${W}s wall for \$N genomes = \$(python3 -c "print(f'{\$W/\$N:.2f}')")s/genome ==="
    )

    export CHUNK=${chunk}   # i/K
    mkdir -p control_genomes   # none in pilot mode
    annotate_genomes ${genomes} ${strains} annotate db ${task.cpus}
    annotate_genomes control_genomes ${control_strains} control_annotate db ${task.cpus}
    """

    stub:
    """
    mkdir -p annotate/${chunk.tokenize('/')[0]}00 control_annotate
    touch annotate/${chunk.tokenize('/')[0]}00/${chunk.tokenize('/')[0]}00.gene_counts.tsv
    """
}
