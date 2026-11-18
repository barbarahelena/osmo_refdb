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

The salt-in negative control, as planned, set 72 salt-in strains (69 of them
archaea) against salt-out strains that are almost all bacteria, and could not
separate osmotic strategy from domain. A post hoc comparison within each
domain, on 621 genomes, is more informative. K⁺ uptake families are enriched
in salt-in lineages in both domains. Salt-in bacteria carry about five fewer
compatible-solute families than equally salt-tolerant salt-out bacteria, and
no detected ectABC, betAB or otsA (otsB is found in 9 of 20). Against
salt-out Clostridia, their own class, the contrast reverses, so it depends on
the reference group. In archaea the compatible-solute half
remains untestable, because the panel finds those families in few archaea of
any kind (test 6b).

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

The plan was written before the full data were collected. It was not
registered with a third party, so we call it prespecified, not pre-registered.
See "Plan provenance" below for what was fixed when. Amendments are listed at
the end of this section.

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
| 6b. the same, within one domain (post hoc) | Does the pattern hold against organisms of the same domain? | as in 6 |

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

**Plan provenance.** The plan lives in this repository's history, with these
timestamps (all 2026-11-04, +0100):

| time | event | evidence |
|---|---|---|
| 11:22 to 12:13 | Pilot: 43 randomly sampled strains harvested, annotated and read out, to test the workflow | file times in `releases/v8/phenotype_validation/` |
| 12:15 | Plan committed: baselines, Sets A and B, blocked cross-validation, the salt-in list and its two predictions, and how each outcome would be read | commit `6ecb9b5` (README and `salt_in_taxa.tsv`) |
| 15:18 | Full BacDive sweep finished (about 1.5 h, with gap-fill reruns) | `harvest*.log` |
| 2026-11-05 | First full analysis | commit `340603d` |

What this does and does not show:

- The plan predates the full data and every full-data result. It postdates
  the 43-strain pilot. Its result was known when the plan was written: rho
  +0.38 between breadth and NaCl max over 41 strains (`analysis.log` in the
  pilot folder), close to the full-data value. The headline association was
  therefore previewed, not predicted blind. The pilot had no salt-in strain
  and no genus eligible for Set B, so the within-genus test, the
  cross-validation and the salt-in control had no preview.
- The commit sat on a local branch until 2026-11-14 and was never pushed
  before the analysis ran. Git commit dates are set by the author's machine.
  There is no independent timestamp, which is why the plan is not described
  as pre-registered anywhere.
- Deviations from the plan are the amendments below. The within-domain
  salt-in comparison (test 6b) was designed after the results were known.

**Salt-in assignment criteria.** A lineage is salt-in when direct
measurements show that it balances external salt with molar intracellular
KCl and holds no significant organic osmolyte pool:

| lineage | rank used | evidence |
|---|---|---|
| Halobacteria | class (Archaea) | the classical salt-in group; reviewed in Oren 2008 (*Saline Systems* 4:2) and Oren 2013 (*Front Microbiol* 4:315) |
| Halanaerobiales | order (Bacteria) | intracellular K⁺, Na⁺ and Cl⁻ at medium strength, no organic solutes, in *Halanaerobium praevalens* and *Halobacteroides halobius* (Oren 1986, *Can J Microbiol* 32:4-9) |
| *Salinibacter*, *Salisaeta* | genus (Bacteria) | *S. ruber*: molar K⁺ and Cl⁻, organic solutes in low concentrations only (Oren et al. 2002, *Extremophiles* 6:491-498); *Salisaeta* assigned as its sister genus |

Everything else is salt-out, including halophiles with "Halo" names that use
compatible solutes (*Halomonas*, *Chromohalobacter*, *Halobacillus*) and the
halophilic methanogens. `salt_in_taxa.tsv` lists the genera one by one, with
"high" confidence where the genus itself was measured and "medium" where it
is assigned from its lineage. The primary analysis uses that list. Applying
the lineage rule directly to BacDive's taxonomy instead moves no strain of the
primary set: the genera the list still misses have no strain with both a
NaCl record and a genome.

Two caveats on the criteria. The salt-in label describes the dominant
strategy, not an absence of compatible-solute genes: 38 of 83 Halobacteria
genomes carry the OtsAB trehalose pathway and 60 carry BCCT-type betaine
transporters (Youssef et al. 2014, *ISME J* 8:636-649). And for the
Halanaerobiales the acidic proteome that goes with salt-in in Halobacteria
was not confirmed by genome analysis (Oren 2013).

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
[0.83, 1.21]. By the reading fixed in the plan, the signal survives control for
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

The prespecified comparison was salt-in strains against all salt-out
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

**6b. Salt-in control within one domain** (post hoc, added 2026-11-14).
The comparison above sets archaea against bacteria. The gene-pattern
prediction needs a genome and a lineage, not a NaCl measurement, so
`build_strategy_control.py` takes every BacDive strain with a genome from
all archaea and from the salt-in bacterial lineages: 302 Halobacteria, 299
other archaea and 20 salt-in bacteria after the 1 Mb filter (510 genomes
beyond the main set; one assembly had no FASTA). Tests use genus means
(Mann-Whitney U) and a genus-cluster bootstrap. Numbers are for v10; v9 gives
the same reading (logs in each release's `strategy_control/`).

| mean families detected per genome | n (genera) | K⁺ uptake (of 5) | compatible solutes (of 19) | random (of 174) |
|---|---|---|---|---|
| Halobacteria (salt-in) | 302 (70) | 1.51 | 1.07 | 41 |
| other archaea | 299 (82) | 0.98 | 1.66 | 33 |
| salt-in bacteria | 20 (13) | 3.30 | 3.70 | 63 |
| salt-out bacteria, NaCl max ≥ 15% | 445 (184) | 2.00 | 8.56 | 72 |
| all salt-out bacteria | 6,520 (1,589) | 1.46 | 5.70 | 74 |

- **K⁺ uptake holds in both domains.** Halobacteria carry 0.53 more K⁺
  uptake families than other archaea [0.28, 0.78] (p = 2×10⁻¹⁰). Salt-in
  bacteria carry 1.30 more than equally salt-tolerant salt-out bacteria
  [0.81, 1.88] (p = 2×10⁻⁶), and 1.04 more than salt-out Clostridia, the class
  most of them belong to. In archaea the gain is no longer separable from
  detectability once the random-pool count and genome size are adjusted for
  (+0.14 [−0.20, +0.49] on v10; +0.36 [+0.05, +0.67] on v9).
- **Compatible solutes hold in bacteria, against the right reference.**
  Salt-in bacteria carry 4.9 fewer compatible-solute families than salt-out
  bacteria that tolerate at least 15% NaCl [−6.0, −3.8] (p = 1×10⁻⁴). No
  salt-in bacterium has a detected ectoine gene, `otsA`, `betA`, `betB` or
  `proP`. They do carry the OpuA transporter (17 of 20) and `betL` (10 of 20),
  so uptake is present and the ectoine, betaine and `otsA` synthesis genes
  are absent; `otsB` is called in 9 of 20 (as in 41% of salt-out bacteria).
  Against salt-out Clostridia,
  which carry few solute families themselves (2.3 on average), salt-in
  strains have more (+1.4), not fewer. Lineage alone is therefore the wrong
  reference for this half; salt tolerance has to be matched.
- **In archaea the solute half stays untestable.** Halobacteria and other
  archaea both carry few of the 19 families (1.07 vs 1.66; difference −0.59
  [−1.29, +0.09]). Ectoine synthesis, `proP`, `proX`, `betA` and `otsA` are
  detected in almost no archaeon of either group, so their absence from
  Halobacteria says nothing about strategy. The one family that separates
  the groups is the OpuA transporter: absent from all 302 Halobacteria,
  present in a quarter of other archaea.
- `otsB` is called in 45% of Halobacteria and 1% of other archaea, with
  `otsA` in none. Youssef et al. (2014) found the OtsAB pathway in 38 of 83
  Halobacteria genomes, so these may be real `otsB` genes whose `otsA`
  partner the panel's bacterial model misses. This has not been checked.

The control is therefore a domain-matched test for bacteria (20 strains, 13
genera) and for K⁺ uptake in archaea. It is not a test of compatible-solute
absence in archaea, and the panel should not be described as validated
there.

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

## Results by release

The Results section above reports v8. The same pipeline was rerun on v9 and
on v10 (v9's models with cutoffs calibrated on a separate split; see
`docs/CHANGELOG.md`). Set A is the same 6,264 strains throughout.

| | v8 | v9 | v10 |
|---|---|---|---|
| breadth rho, Set A | +0.364 [+0.292, +0.437] | +0.373 [+0.299, +0.446] | +0.379 [+0.304, +0.453] |
| breadth rho, adjusted for baselines | +0.508 | +0.510 | +0.513 |
| random panels reaching the osmo panel (of 1,000, any statistic) | 0 | 0 | 0 |
| within-genus partial r | +0.123 | +0.116 | +0.119 |
| median within-genus rho | +0.173 | +0.145 | +0.132 |
| mixed model, % NaCl per SD of breadth | +1.02 [0.83, 1.21] | +0.99 [0.80, 1.18] | +1.01 [0.82, 1.20] |
| ridge CV, genera held out: baselines / + random / + osmo | 0.21 / 0.38 / 0.64 | 0.23 / 0.41 / 0.64 | 0.23 / 0.41 / 0.64 |
| ridge CV, families held out: baselines / + random / + osmo | 0.13 / 0.32 / 0.63 | 0.15 / 0.33 / 0.61 | 0.15 / 0.34 / 0.61 |

Every conclusion holds on all three. The median within-genus rho drifts
down; it is a median over about 315 small genera, and the pooled
within-genus partial r, which has the random-panel null, does not move.

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

The workflow runs it, after the gate (`--gate_passed true`, `docs/hpc.md`), on the frozen release. By
default on the 43 pilot strains (`assets/bacdive/pilot_strains.tsv`), a quick check of the panel; with
`--bacdive_full` on every BacDive strain with a NaCl maximum and a genome, with the within-domain salt-in
control and the NaCl model.

| Step | What it does |
|---|---|
| `BACDIVE_HARVEST` | `--bacdive_full` only: BacDive harvest (`harvest_bacdive.sh`, about 1,820 calls, 1.5 h; batches in `--bacdive_snapshot` are not fetched again), the phenotype table (`build_table.py`), the salt-in control's strains (`build_strategy_control.py`) |
| `BACDIVE_GENOMES` | the genomes from NCBI (gzipped, about 9.4 GB for the full set) |
| `POOL_PFAM`, `POOL_EXTRACT` | the random-panel pool: the 174 families of `assets/random_panel/pool_families.tsv`, their models taken from Pfam-A 39.0, so it is the same for every release and set |
| `BACDIVE_PROTEINS` | per genome, in chunks: genome size and GC, Prodigal proteins (`-p meta`, as osmotool annotate runs it; kept for the next step), the random-pool search (about 5.5 h at 8 jobs for the full set) |
| `BACDIVE_ANNOTATE` | per genome, in the same chunks: osmotool annotate on those proteins, so Prodigal is not rerun |
| `BACDIVE_MATRIX` | the feature matrix (`build_matrix.py`; drops assemblies under 1 Mb) |
| `BACDIVE_ANALYZE` | tests 1-6 (`analyze.py`, `analyze_extended.py`) and, with `--bacdive_full`, the salt-in control, test 6b (`analyze_strategy_control.py`) |
| `NACL_MODEL` | `--bacdive_full` only: the NaCl model (`docs/nacl_model.md`) |

Genomes, proteins and the pool search do not depend on the release: a run for another release reuses
them through `-resume`. The release's results are published to `releases/<release>/phenotype_validation/`
(the annotation, the matrix, the test tables and their printed reports, `*.txt`); the BacDive tables,
harvested batches and the pool's models to `<repo>/<data_dir>`. The scripts each step runs are in its module
(`modules/local/bacdive/<step>/`, `modules/local/pool/`).

Notes:

- Run one thread per job and parallelise across genomes. In the pilot, extra
  threads per job made the HMM step slower (measured with hmmscan).
- The feature matrix has 44 columns, one more than `families.yaml`. The extra
  column is `mrpA_mrpB_fused`, from `build_fusion_refs.py`. It is detected
  in no genome here.
- On macOS arm64 with numpy 2.0, `analyze_extended.py` prints many
  `encountered in matmul` warnings. They come from numpy's Accelerate backend,
  and a plain `a @ b` on random numbers triggers them too. The output tables
  are finite. `grep -v -e Warning -e '^  re'` strips the warnings from the log.
