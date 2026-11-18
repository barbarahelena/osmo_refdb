# Settings of the catalog-index study (README.md). Sourced by every step; override any of them in
# the environment, e.g. CATALOG_GENES=/path/to/GMGC10...fna.gz bash run.sh
STUDY_DIR="${STUDY_DIR:-/scratch-shared/${USER}/catalog_index_study}"   # work and results; visible to compute nodes
ENV_DIR="${ENV_DIR:-${STUDY_DIR}/env}"                                    # 00_setup.sh makes it
RELEASE="${RELEASE:-${STUDY_DIR}/release/v9}"                             # an osmo_refdb release (00_setup.sh downloads v9)
CATALOG_GENES="${CATALOG_GENES:-/projects/0/prjs1784/data/gmgc/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz}"
OSMOTOOL_SRC="${OSMOTOOL_SRC:-/gpfs/work4/0/prjs1784/pipelines/osmotool}"
OSMOTOOL_REF="${OSMOTOOL_REF:-catalog-default}"                           # the osmotool commit installed into the env

N_GENES="${N_GENES:-1000000}"          # genes in the slice (of ~11.9 M)
N_PAIRS="${N_PAIRS:-2000000}"          # read pairs per divergence level
DIVERGENCES="${DIVERGENCES:-0 0.02}"   # per-base divergence of the sample's strain from its catalog gene
ERROR_RATE="${ERROR_RATE:-0.005}"      # sequencing substitution errors per base
SEED="${SEED:-20261118}"
THREADS="${THREADS:-32}"
PARTITION="${PARTITION:-genoa}"
