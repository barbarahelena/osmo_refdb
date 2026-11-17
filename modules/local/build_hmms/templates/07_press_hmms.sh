#!/bin/bash
# ---------------------------------------------------------------------------
# osmo_refdb — 07_press_hmms.sh
# Concatenate all per-family HMMs (with GA cutoffs already set by
# 06_calibrate_cutoffs.py) and hmmpress into a single binary HMM database
# usable with a single hmmscan/hmmsearch call.
#
# Usage: bash 07_press_hmms.sh [hmms_dir] [release_name]
# ---------------------------------------------------------------------------

set -euo pipefail

HMM_DIR="${1:-hmms}"
RELEASE_NAME="${2:-osmo_refdb}"
COMBINED="${HMM_DIR}/${RELEASE_NAME}.hmm"

rm -f "${COMBINED}" "${COMBINED}".h3*
cat "${HMM_DIR}"/*.hmm > "${COMBINED}"
hmmpress "${COMBINED}"

echo "Pressed HMM database: ${COMBINED} (+ .h3f/.h3i/.h3m/.h3p)"
