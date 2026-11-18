#!/bin/bash
# The study's environment (osmotool's tools, bwa-mem3, osmotool itself at OSMOTOOL_REF) and the release.
# Runs on a login node: downloads only.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "${HERE}/config.sh"
source ~/miniforge3/etc/profile.d/conda.sh
mkdir -p "${STUDY_DIR}"

if [ ! -x "${ENV_DIR}/bin/bwa-mem3" ]; then
    mamba create -y -q -p "${ENV_DIR}" -c conda-forge -c bioconda \
        python=3.11 diamond=2.2.8 hmmer=3.4 prodigal orfm bwa-mem3 pandas numpy pyyaml matplotlib scikit-learn pip
fi

# osmotool from one commit, not the working tree (so a checkout elsewhere changes nothing here)
src=$(mktemp -d "${STUDY_DIR}/osmotool_src.XXXX")
git -C "${OSMOTOOL_SRC}" archive "${OSMOTOOL_REF}" | tar -x -C "${src}"
commit=$(git -C "${OSMOTOOL_SRC}" rev-parse --short "${OSMOTOOL_REF}")
SETUPTOOLS_SCM_PRETEND_VERSION="0.8.0+${commit}" "${ENV_DIR}/bin/pip" install -q --no-deps --force-reinstall "${src}"
rm -r "${src}"
echo "${commit}" > "${STUDY_DIR}/osmotool_commit.txt"
"${ENV_DIR}/bin/osmotool" --version

if [ ! -f "${RELEASE}/osmo_refdb.dmnd" ]; then
    "${ENV_DIR}/bin/osmotool" download-db --release "$(basename "${RELEASE}")" --location "$(dirname "${RELEASE}")"
fi
