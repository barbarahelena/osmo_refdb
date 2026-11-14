# Phenotype validation

Does a genome's osmo_refdb profile predict how much salt the organism
tolerates, once genome size and lineage are accounted for?

The `benchmark` step and `real_genome_validation/` test detection: whether
DIAMOND and the HMMs find a family when it is present. This tier tests whether
the families the panel finds relate to a measured phenotype. We use the maximum
NaCl concentration at which each strain grows, curated by DSMZ in BacDive. That
measurement was made independently of our annotations, so the test can fail.

## Summary

Across 6,264 strains, genomes with more of the 44 panel families tolerate more
salt. The Spearman rho is +0.36, and +0.51 after adjusting for genome size,
protein count and GC. The association also holds among strains of the same
genus, at a smaller size: six more panel families (one SD) go with about 1
percentage point higher NaCl tolerance. None of 1,000 random panels of
non-osmotic Pfam families matched the osmo panel on any test. When whole
taxonomic families are held out of model training, the panel keeps its
predictive accuracy and random panels lose most of theirs. Ectoine synthesis
(`ectABC`) and the glycine betaine transporter `betL` carry the most signal
within genera.

The salt-in negative control could test only half of its prediction. A salt-in
haloarchaeal genome carries on average 2.1 of the 5 potassium uptake families,
against 0.9 in salt-out archaea, as predicted. Salt-in genomes carry few of the
19 compatible-solute families (2.0), but salt-out archaea carry few as well
(1.3), so
that half cannot separate osmotic strategy from domain.

## Data

We swept every BacDive ID through the v2 API on 2026-11-04. For each strain we
took the highest NaCl concentration with reported growth ("NaCl max", % w/v)
and the genome assembly BacDive links to that same strain. Matching at strain
level matters because salt tolerance varies within species.

| | strains |
|---|---|
| in BacDive | 102,187 |
| with a NaCl record and a genome accession | 6,643 |
| genome downloaded and annotated | 6,637 |
| excluded as partial assemblies (under 1 Mb) | 3 |
| with a quantitative NaCl max | 6,336 |
| salt-in, by `salt_in_taxa.tsv` | 72 |
| salt-out, the primary analysis set | 6,264 in 1,571 genera |

NaCl max has a median of 5% and ranges from 0 to 55%. We annotated genomes
with `osmotool annotate` v0.7.0 (Prodigal `-p meta`, hmmsearch `--cut_ga`).
**Panel breadth** is the number of the 44 feature-matrix families detected in a
genome: the 43 in `families.yaml` plus the fused `mrpA_mrpB` form. Its mean is
18.7 and its SD 6.2.

We considered the Madin et al. trait compilation and did not use it. Its
salinity field is categorical, covers 684 species, and maps "Mesophilic" (a
temperature trait) into the moderate-halophile class.

## Design and expectations

The plan was written before the full data were collected (commit `6ecb9b5`).
Amendments are listed at the end of this section.

Lineage is the main confounder. Genus alone explains 75% of the variance in
NaCl max, and taxonomic family explains 46%. The panel families are also
distributed by lineage. A raw correlation could therefore arise
because halophilic lineages happen to carry more of these genes, for reasons
unrelated to salt. Genome size
is a second confounder, since larger genomes carry more genes of every kind.
Each test below targets one of these.

| test | question | expected if the panel carries salt biology |
|---|---|---|
| 1. baselines | Does breadth predict NaCl max beyond genome size, protein count and GC? | partial correlation stays positive |
| 2. random panels | Would any 44 gene families do as well? | the osmo panel scores above nearly all random panels |
| 3. within genus | Among strains of one genus, do those with more panel families tolerate more salt? | positive within-genus association |
| 4. blocked cross-validation | Does the panel predict NaCl max in lineages the model has not seen? | accuracy holds when whole genera or families are held out |
| 5. family ranking | Which families carry the signal? | known osmotic genes at the top |
| 6. salt-in control | Do salt-in organisms show their known gene pattern? | K⁺ uptake present, compatible solutes absent |

**Random panels (test 2).** We drew candidate families from Pfam 39.0. We
excluded the 44 panel accessions and 1,126 families with osmotic or cation
terms in their name or description. We screened 6,000 candidates against 43
pilot proteomes and kept the 1,106 detected in at least two. From those we drew
a pool of 174, stratified so that its median prevalence (0.42) matches the
panel's. Each random panel is 44 families from this pool. In the
prevalence-matched version, each osmo family is swapped for one of its five
nearest-prevalence pool families. The p-value is the share of 1,000 random
panels that score at least as high as the osmo panel, so its floor is 1/1001.

**Within genus (test 3).** Set B contains genera with at least 3 strains whose
NaCl max spans at least 4 percentage points: 318 genera and 4,261 strains. We
measured within-genus association three ways:

- the median of per-genus Spearman rho
- the partial correlation after subtracting genus means and adjusting for the
  baselines
- a mixed model, `nacl_max ~ breadth + proteins + genome size + GC + (1|genus)`

The reading was committed in advance. If the within-genus association stays
positive and the top families are osmotic genes, the signal reflects salt
biology. If it collapses toward zero while the across-genus correlation stays
high, the signal was lineage.

**Blocked cross-validation (test 4).** We used 5 folds, with all strains of a
genus (or a family) placed in the same fold. Ridge regression and gradient
boosting were trained on three feature sets:

- baselines alone
- baselines plus a random panel (the median over 20 prevalence-matched draws)
- baselines plus the osmo panel

The score is the out-of-fold Spearman correlation and R².

**Salt-in control (test 6).** Salt-in organisms balance external salt by
accumulating molar KCl in the cytoplasm, and make few or no compatible solutes.
The list covers class Halobacteria, the Haloanaerobiales, and the genera
Salinibacter and Salisaeta. The prediction is that K⁺ uptake families (`trk`, `ktr`) are
common in them and compatible-solute families (`ect`, `bet`, `ots`, `pro`,
`opu`) are rare. Membership was fixed in `salt_in_taxa.tsv` before any panel
output existed, and amended once (see below). Salt-in strains are left out of tests 1 to 5. A secondary model
over all strains adds a strategy term and a strategy × breadth interaction.

**Amendments.**
- 2026-11-04. The plan expected Set B to need extra harvesting of multi-strain
  genera. The full sweep already contained over 300 eligible genera, so we built
  Set B from it directly.
- 2026-11-04. We skipped the planned targeted harvest of salt-in genera,
  because the full sweep already holds every BacDive strain.
- 2026-11-05. We specified the random-panel construction described above. The
  plan said "a size-matched random panel of 43 families". We used 1,000 panels
  of 44 families, matching the feature matrix's 44 columns. This was written
  after the osmo panel's own results were known and before any random-panel
  result existed.
- 2026-11-06. We added 14 Halobacteria genera (25 strains) to
  `salt_in_taxa.tsv`. The list named haloarchaeal genera one by one and had
  missed these, so they were first analysed as salt-out. We found them after
  the full analysis, while checking which archaea the salt-in comparison used.
- 2026-11-06. We excluded assemblies under 1 Mb. Three strains in the sweep are
  linked to partial assemblies (5.3 kb, 68 kb and 0.9 Mb in 1,095 contigs);
  every other assembly is at least 1.35 Mb. Fragments like these carry almost
  no panel families and would count as low-breadth genomes.

Neither 2026-11-06 change altered a conclusion. The largest shift in a
reported number was baselines-only cross-validation with genera held out,
from +0.17 to +0.21. The results below include both changes, and the earlier
outputs are kept in `pre_amendment_2026-11-06/` in the output folder.

## Results

Brackets are 95% genus-cluster bootstrap intervals.

**1. Baselines.** In this dataset, salt-tolerant strains have slightly smaller
genomes and lower GC. Adjusting for these baselines therefore raises the
panel's correlation.

| Spearman rho with NaCl max | |
|---|---|
| total proteins | −0.19 [−0.26, −0.13] |
| genome size | −0.21 [−0.27, −0.15] |
| GC | −0.11 [−0.19, −0.02] |
| panel breadth | +0.36 [+0.29, +0.44] |
| panel breadth, adjusted for all three | +0.51 [+0.44, +0.57] |

**2. Random panels.** No random panel reached the osmo panel on any statistic,
in either null (p ≤ 0.001 throughout). Before adjustment, random panels
correlate negatively with NaCl max, because random gene content tracks genome
size. After adjustment they sit near zero.

| statistic | osmo panel | random median | random 95th percentile |
|---|---|---|---|
| rho, all strains | +0.36 | −0.11 | −0.04 |
| rho, adjusted for baselines | +0.51 | +0.06 | +0.15 |
| within-genus partial r | +0.12 | −0.01 | +0.02 |
| median within-genus rho | +0.17 | 0.00 | +0.02 |

Values are from the prevalence-matched null. The unconstrained null gives the
same picture.

**3. Within genus.** The median within-genus rho is +0.17. Of the 312 genera
with a defined rho, 192 are positive and 105 negative (two-sided sign test
p = 5×10⁻⁷). Among the 26 genera with at least 30 strains, the split is 22 to 4.
The partial correlation after removing genus means and baselines is +0.12. In
the mixed model, one SD of breadth adds 1.02 percentage points of NaCl max
[0.83, 1.21]. By the pre-committed reading, the signal survives control for
lineage, though it is much weaker within genera than across them.

**4. Blocked cross-validation** (ridge regression; out-of-fold Spearman, with
R² in brackets). Gradient boosting ranks the feature sets the same way.

| features | genera held out | families held out |
|---|---|---|
| baselines | +0.21 (0.04) | +0.13 (0.01) |
| baselines + random panel | +0.38 (0.15) | +0.32 (0.06) |
| baselines + osmo panel | +0.64 (0.43) | +0.63 (0.39) |

Random families add some accuracy, since gene content in general carries
lineage information. That gain shrinks when whole families are held out. The
osmo panel's accuracy barely changes, so its signal carries over to families it
was not trained on.

**5. Family ranking** (correlation with NaCl max within genus, Set B).

| family | role | within-genus r |
|---|---|---|
| `ectA`, `ectB`, `ectC` | ectoine synthesis | +0.08 to +0.11 |
| `betL` | glycine betaine uptake | +0.10 |
| `mrpA`, `mrpC`–`mrpG` | Na⁺/H⁺ antiporter complex | +0.06 to +0.08 |
| `ktrA` | K⁺ uptake | +0.06 |
| `nhaA` | Na⁺/H⁺ antiporter | −0.03 |
| `kdpA` | high-affinity K⁺ uptake | −0.06 |

The Mrp subunits are encoded together, so they count as one observation
(`mrpB` is lower, +0.02, and detected in fewer genomes).
`kdpA` is also the most negative family across genera (rho −0.35). Kdp scavenges
K⁺ when it is scarce, which fits low-salt habitats. Bulk K⁺ uptake in
halotolerant organisms runs through Trk and Ktr.

**6. Salt-in control.** The prediction concerns three groups of families. For
each group, the table gives the average number of its families detected per
genome:

- K⁺ uptake: 5 families (`trkA`, `trkH`, `ktrA`, `ktrB`, `ktrD`)
- compatible-solute synthesis and transport: 19 families (`ectABC`, `betABL`,
  `otsAB`, `proPX`, `opuA*`, `opuB*`, `opuC*`)
- the 174 random pool families, as a reference for how many families of any
  kind the HMMs find in a genome

The pre-registered comparison was salt-in strains against all salt-out
strains, and there the prediction holds. A salt-in genome carries 2.1 of the
5 K⁺ uptake families, against 1.6 in salt-out genomes. It carries 2.0 of the
19 compatible-solute families, against 6.1.

Salt-out strains are almost all bacteria, while 69 of the 72 salt-in strains
are archaea, so we also compared against the archaea. This comparison was
added after the results were known.

| mean families detected per genome | n | K⁺ uptake (of 5) | compatible solutes (of 19) | random (of 174) |
|---|---|---|---|---|
| salt-in | 72 | 2.1 | 2.0 | 41 |
| salt-out archaea | 29 | 0.9 | 1.3 | 32 |
| salt-out bacteria | 6,234 | 1.6 | 6.2 | 75 |

An archaeal genome yields about half as many random-family hits as a bacterial
one (32 to 41, against 75). The HMMs find fewer families of every kind in
archaea, so salt-out archaea are the fair reference for salt-in strains:

- Salt-in genomes carry more than twice as many K⁺ uptake families (2.1 vs
  0.9), while their random-family count is only about a quarter higher (41 vs
  32). This half of the prediction holds beyond the domain effect.
- Compatible-solute families are scarce in both groups (2.0 and 1.3 of 19), so
  this half cannot be tested with these data.

The secondary model's strategy × breadth interaction is uninformative:
−0.96 [−4.0, +2.0]. The 72 salt-in strains fall in 36 genera, none of which
also holds salt-out strains, so the model cannot separate strategy from genus.

## Conclusions

The osmo panel carries information about salt tolerance beyond genome size and
lineage. Three results support this:

- The association holds within genera.
- No random set of gene families reproduces it.
- It predicts NaCl max in taxonomic families the model never saw.

The effect is large across lineages and modest within a genus. Within a genus,
breadth explains about 1.5% of the remaining variance (r = 0.12). Profiles can
therefore rank genomes by likely salt tolerance, but they will not pin down a
single strain's NaCl max: out-of-fold R² is about 0.4. This study used
assemblies only. How well the result carries over to read-based community
profiles is untested.

The informative families are ectoine synthesis, BetL and the Mrp antiporter,
with `kdpA` pointing the other way. A score that weights families, or reports
these groups separately, may outperform raw breadth. We have not tested this.

In archaeal genomes the panel detects fewer families of every kind, and its
compatible-solute families say little about archaea. Breadth is likely to
read low for archaea-rich samples such as hypersaline brines, and archaeal
profiles should be read separately from bacterial ones.

## Limitations

- Genus names follow phylogeny unevenly. Blocking on 16S identity or tree
  distance would be stricter, and BacDive records carry 16S accessions for
  this.
- NaCl max comes from many labs, media and test ranges. It is also censored at
  the highest concentration each lab tested.
- All random panels come from one pool of 174 families and overlap. Their null
  is narrower than independent draws from all of Pfam would give.
- 163 NaCl records in molar or g/L units are not yet parsed.
- Only 29 salt-out archaea are available as a reference for the salt-in
  comparison, from 21 genera.

## Running it

Scripts live here. Output goes to `releases/<release>/phenotype_validation_full/`
(git-ignored). Use osmotool v0.7.0 or later: older pins lack part of the
43-family list, and `annotate` under-reports without it.

```bash
PV=releases/v8/phenotype_validation_full
RP="$PV/random_panel"

# 1. BacDive harvest (~1,820 calls, ~1.5 h). Cached batches are skipped, so a
#    rerun only fetches what failed. Keep sweep and sample runs in separate dirs.
bash  pipeline/phenotype_validation/harvest_bacdive.sh "$PV" sweep

# 2. phenotype table, one row per strain with NaCl max and a genome accession
python pipeline/phenotype_validation/build_table.py --dir "$PV"

# 3. genomes from NCBI (kept gzipped, ~9.4 GB)
PAR=8 bash pipeline/phenotype_validation/fetch_genomes.sh "$PV"

# 4. osmotool annotate in the pinned image, one thread per job
docker run --rm -v "$PWD":/work -w /work --entrypoint bash osmo_refdb:latest \
    pipeline/phenotype_validation/run_annotate.sh "$PV" releases/v8 10 1

# 5. feature matrix (drops assemblies under 1 Mb) and first readout
python pipeline/phenotype_validation/build_matrix.py --dir "$PV"
python pipeline/phenotype_validation/analyze.py      --dir "$PV"

# 6. random-family pool from Pfam-A.hmm.gz, screened on the 43 pilot proteomes
python pipeline/phenotype_validation/build_random_pool.py candidates \
    --pfam "$RP/Pfam-A.hmm.gz" --families families.yaml --out-dir "$RP" --n 6000
hmmsearch --cut_ga --cpu 8 --noali --tblout "$RP/screen.tbl" -o /dev/null \
    "$RP/candidates.hmm" "$RP/screen_proteomes.faa"
python pipeline/phenotype_validation/build_random_pool.py pool \
    --pfam "$RP/Pfam-A.hmm.gz" --out-dir "$RP" --screen "$RP/screen.tbl" \
    --osmo-matrix releases/v8/phenotype_validation/feature_matrix.tsv --n-genomes 43

# 7. genome size, GC and pool search per genome (~5.5 h at 8 jobs), then
#    tests 1-6 (needs numpy, pandas, scipy, scikit-learn, statsmodels)
docker run --rm -v "$PWD":/work -w /work --entrypoint bash osmo_refdb:latest \
    pipeline/phenotype_validation/run_random_panel.sh "$PV" "$RP/pool.hmm" 8
python pipeline/phenotype_validation/analyze_extended.py --dir "$PV"
```

Notes:

- Run one thread per job and parallelise across genomes. In the pilot, extra
  threads per job made the HMM step slower (measured with hmmscan).
- The feature matrix has 44 columns, one more than `families.yaml`. The extra
  column is `mrpA_mrpB_fused`, from `08d_build_fusion_refs.py`. It is detected
  in no genome here.
- On macOS arm64 with numpy 2.0, `analyze_extended.py` prints many
  `encountered in matmul` warnings. They come from numpy's Accelerate backend,
  and a plain `a @ b` on random numbers triggers them too. The output tables
  are finite. `grep -v -e Warning -e '^  re'` strips the warnings from the log.
