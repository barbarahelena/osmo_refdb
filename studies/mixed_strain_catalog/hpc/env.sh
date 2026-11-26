# Settings for running the study on a SLURM cluster with Apptainer. Sourced by
# prepare.sh and the two jobs; edit the first block (or export the variables first).
#
# --- edit
export RELEASE=${RELEASE:-/path/to/releases/v11}   # the release whose shipped index is compared
export CATALOG_GENES=${CATALOG_GENES:-/path/to/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz}
                                                    # the build's --catalog_genes; prepare.sh fetches it there if missing
export FULL_INDEX=${FULL_INDEX:-/path/to/gmgc_gut_full_index}
                                                    # built once by full.sbatch, reused if present (profile --catalog_index)
export SIF_DIR=${SIF_DIR:-$HOME/containers}         # where prepare.sh puts the two images
APPTAINER=${APPTAINER:-apptainer}
# --- derived
for v in RELEASE CATALOG_GENES FULL_INDEX; do
    case ${!v} in /path/to/*) echo "hpc/env.sh: set $v (still ${!v})" >&2; return 1 2>/dev/null || exit 1 ;; esac
done
STUDY=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=$(cd "$STUDY/../.." && pwd)
export STUDY
export CATALOG=${CATALOG:-$RELEASE/catalogs/gmgc_gut}
export THREADS=${SLURM_CPUS_PER_TASK:-${THREADS:-8}}
OSMOTOOL_SIF=$SIF_DIR/osmotool_0.11.1.sif               # osmotool 0.11.1 + bwa-mem3 0.18.0, as run locally
WGSIM_SIF=$SIF_DIR/wgsim_1.0--h577a1d6_10.sif            # the pipeline's WGSIM image
mkdir -p "$STUDY/work" "$(dirname "$CATALOG_GENES")" "$(dirname "$FULL_INDEX")"
BINDS=$(for d in "$ROOT" "$RELEASE" "$(dirname "$CATALOG_GENES")" "$(dirname "$FULL_INDEX")"; do
            readlink -f "$d"; done | sort -u | paste -sd, -)
W=$STUDY/work
export OSMOTOOL_RUN="$APPTAINER exec --bind $BINDS --pwd $W $OSMOTOOL_SIF bash -c"
export BWA_RUN="$OSMOTOOL_RUN"
export WGSIM_RUN="$APPTAINER exec --bind $BINDS --pwd $W $WGSIM_SIF bash -c"
export PYTHON_RUN="$APPTAINER exec --bind $BINDS --pwd $W $OSMOTOOL_SIF python3"
