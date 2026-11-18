# Tests

Two pipeline tests (`main.nf.test`, nf-test) on the `test` profile's inputs (`assets/test/`):

| test | what it checks | how long | run with |
|---|---|---|---|
| `stub` | the wiring: every step's stub runs, through the gate, and the release has its layout | about 3 min, offline | `NXF_OFFLINE=true nf-test test --tag stub` |
| `real` | a whole release built for real on small inputs, then `check_performance.py` | about an hour, needs the network | `nf-test test --tag real --profile test,snellius` |

The default profile is `test,singularity`: the test inputs, containers on the machine nf-test runs
on, no SLURM. Even a stub runs in each step's container, because its version outputs read the
container's tools. On a cluster, run nf-test itself inside a SLURM job, like any Nextflow run.
`NFT_WORKDIR` moves nf-test's work directory (default `.nf-test/`) to scratch;
`NXF_SINGULARITY_CACHEDIR` keeps the images between runs.

## check_performance.py

It scores a release, the test one or any other (`--release-dir releases/<name>
--families families.yaml --markers markers.yaml`). Standard library only.

- **Must hold exactly:**
  - the release's parts are there;
  - every family has an HMM with a gathering cutoff, a reporting unit and a benchmark row;
  - every marker is in `markers/`;
  - the per-read rebuilds agree with osmotool's own counts, in `results/metrics/verification.tsv`
    and `results/catalog_benchmark/verification.tsv`;
  - the catalog benchmark's median f1_all with every sequence catalogued (A) is at least that with
    none of the scored ones (B).
- **Must not drop by more than `--tolerance` (0.05) below `baseline.tsv`:**
  - the read benchmark: f1_all, precision_all and recall per reporting unit, for diamond_cascade and hmm;
  - the annotate benchmark: f1, precision and recall per family;
  - the catalog benchmark: per unit, A and B;
  - the real-genome reads: recall and precision per target family, without and with the catalog.
- **Every run** writes its metrics to `performance.tsv`.

The references come from UniProt as it is that day, so a release built from another UniProt release
than the baseline's can differ for that reason alone. The baseline records its UniProt release, and
across releases a drop is reported as a warning rather than a failure (`--strict` makes it fail).

### The baseline

`baseline.tsv` is written by the first real run, or with `--write-baseline`, from that run's metrics.
Review it before committing it: it is what later builds are held to. Rewrite it when a change is meant
to move the numbers, and say so in `CHANGELOG.md`.

The reference-system check (`real_genome_validation/reference_systems_check.txt`) is printed but not
scored: its expectations are for the full panel, not the test's five families.
