// The canonical build of a release: run_pipeline.sh all from a committed checkout,
// inside the pinned image, then PROVENANCE.json. From the snapshot at reuse_refs_from
// everything after the fetch (split, alignment of changed families, HMMs, calibration,
// DIAMOND, benchmark) is rerun. No automatic retry: a failed build leaves a partial
// releases/<name>; remove it before running again.
process BUILD_RELEASE {
    // calls: the whole toolchain AND osmotool: mafft, trimal, hmmer, orfm, wgsim, diamond, cd-hit, prodigal, bwa-mem3, bwa, insilicoseq, python + requests, pyyaml, numpy, pandas, scikit-learn, scipy, statsmodels, matplotlib, xgboost, biopython; curl, git, procps
    // image: see hpc/containers/README.md: osmotool + the extra build tools in one image
    container 'todo.invalid/build_release'
    label 'process_high'
    input:
    tuple val(commit), val(dirty), val(diff)

    output:
    val(params.release)

    script:
    """
    cd ${params.repo}
    if [ -e releases/${params.release} ]; then
        echo "ERROR: releases/${params.release} exists; remove it (a failed earlier build?) or use a new release name" >&2
        exit 1
    fi
    if [ "${dirty}" != "0" ]; then
        echo "WARNING: ${dirty} tracked files differ from commit ${commit}; PROVENANCE.json will list them" >&2
    fi
    echo "commit ${commit}  uncommitted-diff ${diff}  node \$(hostname)  arch \$(uname -m)  cores with AVX2: \$(grep -c avx2 /proc/cpuinfo || true)"
    THREADS=${task.cpus} JOBS=${task.cpus} REUSE_REFS_FROM=${params.reuse_refs_from} SPLIT_MODE=${params.split_mode} \\
        bash run_pipeline.sh all ${params.release}
    python pipeline/write_provenance.py --release releases/${params.release}
    test -s releases/${params.release}/hmms/osmo_refdb.hmm
    """

    stub:
    """
    echo stub
    """
}
