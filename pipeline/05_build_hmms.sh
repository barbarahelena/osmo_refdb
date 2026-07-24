#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 05_build_hmms.sh
# hmmbuild from each family's trimmed TRAIN alignment, then score the
# positive test set and hard-negative set against it (used by
# 06_calibrate_cutoffs.py). Assumes hmmer is on PATH.
#
# NOTE: hmmbuild is given -n "${FAMILY}" explicitly so the resulting HMM's
# model name is the family code (e.g. "ectC"), not the alignment file's
# internal MSA name (e.g. "ectC.trimmed") — otherwise downstream family
# matching in compute_metrics.py silently fails (all HMM calls attributed
# to the wrong "family").
#
# families.yaml `pfam_model: PFxxxxx` (added 2026-07-23, Task 1 gap-filling
# genes): for a family whose true biological target IS essentially a single
# curated Pfam family (checked per-gene, not assumed panel-wide -- see
# families.yaml comments for which of the new genes qualify and which
# don't), fetch that family's actual, curated Pfam-A HMM -- built from
# Pfam's own hand-curated seed alignment, carrying its own GA/TC/NC
# gathering-threshold lines -- instead of hmmbuild-ing a fresh model from
# our own MAFFT/trimAl alignment. Only the model NAME line is rewritten (to
# our family code, for the same reason noted above); the GA cutoff Pfam
# shipped is left untouched, not recomputed -- see
# 06_calibrate_cutoffs.py for how it's then used (read, not overwritten).
# Confirmed working via the allowed www.ebi.ac.uk host:
#   curl ".../interpro/api/entry/pfam/<PFxxxxx>/?annotation=hmm" returns a
#   gzipped, valid HMMER3 .hmm file including its GA line -- NOT the same
#   as the InterPro entry metadata endpoint's `entry_annotations.hmm`
#   counter, which reports 0 regardless (that counter tracks something
#   else; don't use it to decide whether a model is fetchable).
# This family still gets a normal positive/negative UniProt fetch (01) for
# DIAMOND references and for scoring below -- only the HMM's SOURCE changes.
#
# Usage: bash 05_build_hmms.sh [refs_dir] [alignments_dir] [hmms_dir] [families_yaml]
# ---------------------------------------------------------------------------

set -euo pipefail

REFS_DIR="${1:-refs}"
ALN_DIR="${2:-alignments}"
HMM_DIR="${3:-hmms}"
FAMILIES_FILE="${4:-families.yaml}"
THREADS="${THREADS:-4}"

mkdir -p "${HMM_DIR}/scores"

pfam_model_for_family() {
    python3 -c "
import sys, yaml
with open(sys.argv[1]) as fh:
    data = yaml.safe_load(fh)
for fam in data['families']:
    if fam['name'] == sys.argv[2]:
        print(fam.get('pfam_model', ''))
        break
" "${FAMILIES_FILE}" "$1"
}

fetch_pfam_hmm() {
    # requests, not curl -- curl isn't installed in the osmo_refdb Docker
    # image (only procps/git are apt-get'd; requests is already a hard
    # dependency of 01_fetch_refs.py, so no new package requirement here).
    local ACCESSION="$1" FAMILY="$2" HMM_FILE="$3"
    python3 -c "
import gzip, io, re, sys
import requests
accession, family, hmm_file = sys.argv[1], sys.argv[2], sys.argv[3]
r = requests.get(f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/', params={'annotation': 'hmm'}, timeout=30)
r.raise_for_status()
text = gzip.decompress(r.content).decode()
text = re.sub(r'^NAME  .*$', f'NAME  {family}', text, count=1, flags=re.MULTILINE)
# hmmpress's SSI index requires unique NAME *and* ACC across the combined
# db -- Pfam's own ACC line (e.g. 'ACC   PF02386.23') is untouched by the
# NAME rewrite above, so two families adopting the SAME pfam_model
# accession (e.g. trkH/ktrB/ktrD all -> PF02386) collide at press time
# unless ACC is also made unique per family. Rewritten to the family name
# too, same as NAME -- families.yaml's own name field is already
# required-unique, so this can't collide.
if re.search(r'^ACC   ', text, flags=re.MULTILINE):
    text = re.sub(r'^ACC   .*$', f'ACC   {family}', text, count=1, flags=re.MULTILINE)
with open(hmm_file, 'w') as fh:
    fh.write(text)
" "${ACCESSION}" "${FAMILY}" "${HMM_FILE}"
}

# Build the list of families to process from families.yaml directly (not
# just alignments/*.trimmed.fasta) so a pfam_model family with too few
# positives to align (04_align_trim.sh's ">=2 sequences" floor) still gets
# processed -- it doesn't need OUR alignment at all.
ALL_FAMILIES=$(python3 -c "
import sys, yaml
with open(sys.argv[1]) as fh:
    data = yaml.safe_load(fh)
print('\n'.join(fam['name'] for fam in data['families']))
" "${FAMILIES_FILE}")

for FAMILY in ${ALL_FAMILIES}; do
    TRIMMED_ALN="${ALN_DIR}/${FAMILY}.trimmed.fasta"
    HMM_FILE="${HMM_DIR}/${FAMILY}.hmm"
    POS_SCORES="${HMM_DIR}/scores/${FAMILY}.positive.tblout"
    NEG_SCORES="${HMM_DIR}/scores/${FAMILY}.negative.tblout"

    PFAM_MODEL=$(pfam_model_for_family "${FAMILY}")
    if [ -n "${PFAM_MODEL}" ]; then
        echo "[${FAMILY}] fetching curated Pfam model ${PFAM_MODEL} (pfam_model in families.yaml)..."
        fetch_pfam_hmm "${PFAM_MODEL}" "${FAMILY}" "${HMM_FILE}"
    else
        [ -s "${TRIMMED_ALN}" ] || { echo "[${FAMILY}] SKIP: no trimmed alignment (run 04 first)"; continue; }
        echo "[${FAMILY}] hmmbuild..."
        hmmbuild --amino -n "${FAMILY}" --cpu "${THREADS}" "${HMM_FILE}" "${TRIMMED_ALN}" > /dev/null
    fi

    # Score against the held-out positive TEST set (not train) as a sanity
    # check of true-positive recall on unseen data.
    if [ -s "${REFS_DIR}/${FAMILY}.positive.test.faa" ]; then
        hmmsearch --cpu "${THREADS}" --tblout "${POS_SCORES}" --noali \
            "${HMM_FILE}" "${REFS_DIR}/${FAMILY}.positive.test.faa" > /dev/null
    elif [ -s "${REFS_DIR}/${FAMILY}.positive.faa" ]; then
        hmmsearch --cpu "${THREADS}" --tblout "${POS_SCORES}" --noali \
            "${HMM_FILE}" "${REFS_DIR}/${FAMILY}.positive.faa" > /dev/null
    fi

    if [ -s "${REFS_DIR}/${FAMILY}.negative.faa" ]; then
        hmmsearch --cpu "${THREADS}" --tblout "${NEG_SCORES}" --noali \
            "${HMM_FILE}" "${REFS_DIR}/${FAMILY}.negative.faa" > /dev/null
    else
        : > "${NEG_SCORES}"   # empty; 06_calibrate_cutoffs.py handles this gracefully
    fi

    echo "[${FAMILY}] done -> ${HMM_FILE}"
done

echo ""
echo "hmmbuild + scoring complete. Next: python 06_calibrate_cutoffs.py"
