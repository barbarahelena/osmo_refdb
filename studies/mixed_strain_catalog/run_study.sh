#!/usr/bin/env bash
# Mixed-strain catalog study: one simulated sample pooling reads from genomes whose
# genes are in the catalog, partly in it, or absent from it (genomes.tsv), scored per
# group for DIAMOND alone and for the catalog path under several settle rules.
#
#   bash run_study.sh [step ...]     steps: genomes annotate reads truth diamond map score verify
#                                    (no argument: all, skipping steps whose output exists)
#
# Environment (defaults are this machine's; on the HPC point them at Apptainer images):
#   RELEASE       osmo_refdb release directory (default ../../releases/v11)
#   CATALOG       the release's catalog directory: catalog.annotation.tsv, the hmmsearch
#                 domtblout and its own (small) index (default $RELEASE/catalogs/gmgc_gut)
#   CATALOG_INDEX directory of the index reads are mapped to (default $CATALOG). For the
#                 full catalog: the directory osmotool index-catalog wrote (catalog.fna.*,
#                 catalog_index.json); genes it holds without an annotation row are background.
#   CATALOG_TAG   name for this index's outputs, work/<CATALOG_TAG>/ (default: the index
#                 directory's name), so runs with different indexes sit side by side
#   COVERAGE      read coverage per genome (default 2)
#   WGSIM_ERR     wgsim base error rate (default 0.02, wgsim's default, as the release's
#                 real-genome validation)
#   THREADS       default 8
#   OSMOTOOL_RUN  how to run osmotool 0.11.1 (default: its Docker image, linux/amd64)
#   NATIVE_RUN    how to run wgsim / bwa-mem3 (default: the osmo_refdb:bwa Docker image)
set -euo pipefail
cd "$(dirname "$0")"
HERE=$PWD
ROOT=$(cd ../.. && pwd)
RELEASE=${RELEASE:-$ROOT/releases/v11}
CATALOG=${CATALOG:-$RELEASE/catalogs/gmgc_gut}
CATALOG_INDEX=${CATALOG_INDEX:-$CATALOG}
CATALOG_TAG=${CATALOG_TAG:-$(basename "$CATALOG_INDEX")}
COVERAGE=${COVERAGE:-2}
WGSIM_ERR=${WGSIM_ERR:-0.02}
THREADS=${THREADS:-8}
W=$HERE/work
MOUNTS="-v $ROOT:$ROOT -v $CATALOG:$CATALOG -v $CATALOG_INDEX:$CATALOG_INDEX"
OSMOTOOL_RUN=${OSMOTOOL_RUN:-"docker run --rm --platform linux/amd64 $MOUNTS -w $W --entrypoint bash docker.io/barbarahelena/osmotool:0.11.1 -c"}
NATIVE_RUN=${NATIVE_RUN:-"docker run --rm $MOUNTS -w $W --entrypoint bash osmo_refdb:bwa -c"}
STEPS=${*:-genomes annotate reads truth diamond map score verify}
mkdir -p "$W"/{genomes,annotate,reads}

for step in $STEPS; do
case $step in
genomes)   # genome FASTAs, contig names prefixed "<genome>__" so every read names its genome
    tail -n +2 genomes.tsv | while IFS=$'\t' read -r g _grp _sp src; do
        out=$W/genomes/$g.fna; [ -s "$out" ] && continue
        case $src in
            bacdive:*) gunzip -c "$ROOT"/releases/v8/phenotype_validation_full/genomes/"${src#bacdive:}"_*.fna.gz ;;
            *)         cat "$ROOT/$src" ;;
        esac | sed "s/^>/>${g}__/" > "$out"
    done ;;
annotate)  # truth genes: osmotool annotate (hmm) with the release, on each genome
    $OSMOTOOL_RUN "cd annotate && for f in ../genomes/*.fna; do g=\$(basename \$f .fna);
        [ -s \$g.gene_coordinates.tsv ] && continue;
        osmotool annotate $RELEASE \$f --out_prefix \$g --threads $THREADS > \$g.stdout 2>&1; done" ;;
reads)     # wgsim per genome (seed = row number), COVERAGE x, pooled into one sample
    [ -s "$W/reads/sample_R1.fq" ] && continue
    i=0; : > "$W/reads/pairs.tsv"
    tail -n +2 genomes.tsv | while IFS=$'\t' read -r g _rest; do
        i=$((i + 1)); bp=$(grep -v '>' "$W/genomes/$g.fna" | tr -d '\n' | wc -c)
        n=$(( bp * COVERAGE / 300 )); printf '%s\t%s\t%s\n' "$g" "$i" "$n" >> "$W/reads/pairs.tsv"
    done
    $NATIVE_RUN "cd reads && while IFS=\$'\t' read g s n; do
        wgsim -S \$s -N \$n -1 150 -2 150 -e $WGSIM_ERR ../genomes/\$g.fna \$g.R1.fq \$g.R2.fq > /dev/null; done < pairs.tsv;
        cat \$(cut -f1 pairs.tsv | sed 's/\$/.R1.fq/') > sample_R1.fq; cat \$(cut -f1 pairs.tsv | sed 's/\$/.R2.fq/') > sample_R2.fq" ;;
truth)     # read truth: overlap >= 20 bp with an annotated panel gene (real_read_truth.py's rule)
    python3 -I "$HERE/truth.py" "$W" ;;
diamond)   # DIAMOND for every read: osmotool profile without a catalog, the validation's options
    [ -s "$W/diamond/sample.blastx.tsv" ] && continue
    mkdir -p "$W/diamond"
    $OSMOTOOL_RUN "cd diamond && osmotool profile $RELEASE -1 ../reads/sample_R1.fq -2 ../reads/sample_R2.fq \
        --out_prefix sample --threads $THREADS --min_identity 0.40 --min_query_cover 0.50 \
        --keep_aln --no_markers --translated_only > sample.stdout 2>&1" ;;
map)       # the catalog's placements (primary + XA), as osmotool profile --catalog maps them
    mkdir -p "$W/$CATALOG_TAG"; [ -s "$W/$CATALOG_TAG/sample.sam" ] && continue
    $NATIVE_RUN "bwa-mem3 mem -v 1 -t $THREADS -h 50 $CATALOG_INDEX/catalog.fna reads/sample_R1.fq reads/sample_R2.fq \
        > $CATALOG_TAG/sample.sam 2> $CATALOG_TAG/bwa.log" ;;
score)     # settle rules x groups
    python3 -I "$HERE/score.py" --work "$W" --release "$RELEASE" --catalog "$CATALOG" --tag "$CATALOG_TAG" \
        --compare "$ROOT/modules/local/validation/readrecall/templates/real_genome_validation/compare_real_recall.py" ;;
verify)    # osmotool profile --catalog itself on the sample: its settled mates must equal score.py's current rule
    mkdir -p "$W/$CATALOG_TAG/osmotool"
    idx=""; [ "$CATALOG_INDEX" != "$CATALOG" ] && idx="--catalog_index $CATALOG_INDEX"
    [ -s "$W/$CATALOG_TAG/osmotool/sample.catalog_hits.tsv" ] || \
    $OSMOTOOL_RUN "cd $CATALOG_TAG/osmotool && osmotool profile $RELEASE -1 ../../reads/sample_R1.fq -2 ../../reads/sample_R2.fq \
        --out_prefix sample --threads $THREADS --min_identity 0.40 --min_query_cover 0.50 --keep_aln --no_markers \
        --catalog_dir $CATALOG $idx > sample.stdout 2>&1"
    python3 -I - "$W/$CATALOG_TAG" <<'EOF'
import csv, sys
d = sys.argv[1]
key = lambda r, rid: (rid, r["decision"], r["label"] if r["decision"] == "call" else "")
mine = {key(r, r["read_id"]) for r in csv.DictReader(open(f"{d}/settled_current.tsv"), delimiter="\t")}
theirs = {key(r, f"{r['read']}/{r['mate']}") for r in csv.DictReader(open(f"{d}/osmotool/sample.catalog_hits.tsv"), delimiter="\t")}
print(f"verify: {len(theirs)} mates settled by osmotool, {len(mine ^ theirs)} rows differ from score.py's current rule")
sys.exit(1 if mine ^ theirs else 0)
EOF
    ;;
*) echo "unknown step $step" >&2; exit 1 ;;
esac
done
