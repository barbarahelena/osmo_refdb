# Catalog-based profiling: benchmark on the v10 test reads

**Conclusion.** When a sample's strains are in the catalog, catalog-first
profiling (`osmotool profile --catalog`) lifts the median f1_all of the 27
reported units from 0.770 to 0.892. Precision rises from 0.890 to 0.989 and
false positives on unrelated reads drop more than tenfold. When the strains
are not in the catalog, nothing is lost: the translated fallback keeps the
v10 numbers (0.800 in this simulation, against 0.770). The exception is
betL, which drops from 0.87 to 0.46 because its full-length HMM cannot
separate betL from betT/caiT. Five excluded units become reportable from
catalog reads alone. The opu ATPases, opu substrate-binding proteins and the
trk/ktr channel do not.

Script: `modules/local/catalog/package/templates/benchmark_catalog.py`. Results:
`releases/v10/results/catalog_benchmark/`. Log:
`releases/v10_catalog_benchmark.log`. osmotool: branch `catalog-profile`
(`df87119`), bwa-mem3 0.17.0.

## Design

The catalog's value depends on whether a sample's strains are in it, so two
bounds are scored on the same 854,484 test reads as steps 9-11:

| scenario | catalog contains | genes | what it bounds |
|---|---|---|---|
| A | every sequence of the release: train, calib and test positives, both negative splits, the random background | 148,721 | upper: every test read's own gene is catalogued |
| B | only what the benchmark does not score: positive.train, positive.calib, negative.train | 105,996 | lower: no test read's own gene is catalogued |

- Both are indexed in full (competitive mapping, as Cayman). A is also
  scored with the small index a release ships (A_hits: genes with an HMM
  hit only).
- Catalog genes are the proteins back-translated exactly as
  `simulate_reads.py` builds its constructs, so a test read matches its
  own gene base for base (checked on 95 constructs).
- Catalog thresholds are Cayman's: >= 97% identity over >= 45 bp aligned.
  Catalogs are annotated with v10's HMMs and gathering cutoffs.
- The fallback reads get v10's DIAMOND + cascade calls, and units are scored
  with v10's report levels, using the same metrics as step 11
  (`summary_by_report_unit.tsv`).

**Per-read rebuild, verified.** Each step decides each read on its own, so
the hybrid call for a read is the catalog's decision when the catalog
settles it, and otherwise the read's v10 translated call. Full `profile
--catalog` runs on three samples (opuAC.negative, proX.positive,
trkH.negative) reproduce the rebuilt counts exactly: 39/39 families each.

## Results

Median over the 27 reported units:

| method | f1_all | precision_all | recall |
|---|---|---|---|
| no catalog (v10) | 0.770 | 0.890 | 0.674 |
| **A, full index** | **0.892** | **0.989** | **0.852** |
| A, release index (hits) | 0.892 | 0.989 | 0.852 |
| A, catalog only (no fallback) | 0.890 | 0.989 | 0.846 |
| B, full index | 0.800 | 0.915 | 0.711 |

Mixture: samples where a fraction f of the reads comes from catalogued
strains (counts mixed linearly between A and B). Median f1_all over the
reported units:

| f | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| median f1_all | 0.800 | 0.831 | 0.849 | 0.867 | 0.892 |

False positives on the random background (unrelated bacterial genes), per
million background reads:

| no catalog | A full | A hits | A catalog only | B full |
|---|---|---|---|---|
| 2083 | 170 | 270 | 170 | 1663 |

What the catalog settled (test reads):

| | A full | A hits | B full |
|---|---|---|---|
| positive reads settled | 106,333 | 106,220 | 28,743 |
| of which called, right family | 98,927 | 98,927 | 27,580 |
| of which called, other family | 492 | 492 | 104 |
| hard-negative reads settled | 406,705 | 282,039 | 74,787 |
| of which called as a panel family | 31,175 | 31,183 | 8,261 |
| background reads settled | 55,081 | 2,252 | 3,771 |

- Settled positive reads are almost always the right family (99.5%).
- About 7.7% of settled hard-negative reads are called as a panel family.
  That is the full-length HMM's hard-negative call rate (annotate: 7.3%).
  Part of it is likely unlabelled real orthologs in the negative pools
  (documented for trkH and otsA in `families.yaml`).

### Per reporting unit (f1_all unless noted)

| unit | level | no catalog | A (full index) | A precision_all | A recall | A, catalog only | A (release hits index) | B (not catalogued) |
|---|---|---|---|---|---|---|---|---|
| ectA | gene | 0.805 | 0.919 | 1.0 | 0.85 | 0.915 | 0.919 | 0.816 |
| ectB | gene | 0.704 | 0.976 | 0.98 | 0.973 | 0.975 | 0.976 | 0.772 |
| ectC | gene | 0.803 | 0.892 | 1.0 | 0.805 | 0.89 | 0.892 | 0.823 |
| teaA | gene | 0.76 | 0.963 | 1.0 | 0.929 | 0.959 | 0.963 | 0.804 |
| betL | gene | 0.868 | 0.459 | 0.465 | 0.454 | 0.454 | 0.459 | 0.776 |
| kdpA | gene | 0.977 | 0.989 | 1.0 | 0.978 | 0.985 | 0.989 | 0.98 |
| nhaA | gene | 0.949 | 0.984 | 1.0 | 0.968 | 0.98 | 0.984 | 0.951 |
| proX | gene | 0.55 | 0.69 | 0.607 | 0.8 | 0.689 | 0.69 | 0.588 |
| opuAC | none | 0.034 | 0.022 | 0.012 | 0.167 | 0.021 | 0.022 | 0.041 |
| opuBC | none | 0.014 | 0 | 0.0 | 0.0 | 0 | 0 | 0.016 |
| opuCC | none | 0.044 | 0.086 | 0.048 | 0.408 | 0.088 | 0.086 | 0.052 |
| proP | gene | 0.494 | 0.682 | 0.673 | 0.691 | 0.68 | 0.688 | 0.505 |
| otsA | none | 0.554 | 0.692 | 0.751 | 0.642 | 0.686 | 0.692 | 0.555 |
| otsB | gene | 0.795 | 0.839 | 0.826 | 0.852 | 0.836 | 0.839 | 0.8 |
| mscL | gene | 0.767 | 0.885 | 0.971 | 0.813 | 0.881 | 0.885 | 0.812 |
| mscS | none | 0.414 | 0.644 | 0.539 | 0.798 | 0.644 | 0.644 | 0.439 |
| galE | gene | 0.858 | 0.932 | 0.98 | 0.889 | 0.929 | 0.931 | 0.869 |
| mazG | gene | 0.888 | 0.953 | 0.997 | 0.913 | 0.95 | 0.953 | 0.895 |
| murB | none | 0.927 | 0.976 | 1.0 | 0.954 | 0.972 | 0.976 | 0.931 |
| opuAA | none | 0.2 | 0.359 | 0.224 | 0.913 | 0.359 | 0.357 | 0.25 |
| opuBA | none | 0.066 | 0.369 | 0.242 | 0.778 | 0.373 | 0.369 | 0.084 |
| opuCA | none | 0.132 | 0.32 | 0.285 | 0.365 | 0.321 | 0.236 | 0.162 |
| opu_permease | group | 0.509 | 0.551 | 0.408 | 0.848 | 0.549 | 0.551 | 0.514 |
| mrpA | none | 0.426 | 0.845 | 0.755 | 0.959 | 0.845 | 0.845 | 0.558 |
| mrpD | none | 0.264 | 0.387 | 0.255 | 0.803 | 0.389 | 0.387 | 0.294 |
| mrpB | gene | 0.716 | 0.838 | 0.868 | 0.81 | 0.838 | 0.838 | 0.759 |
| mrpC | none | 0.503 | 0.572 | 0.472 | 0.724 | 0.568 | 0.571 | 0.517 |
| mrpE | gene | 0.691 | 0.929 | 1.0 | 0.867 | 0.929 | 0.929 | 0.696 |
| mrpF | gene | 0.704 | 0.882 | 1.0 | 0.789 | 0.882 | 0.882 | 0.718 |
| mrpG | gene | 0.68 | 0.903 | 1.0 | 0.824 | 0.902 | 0.903 | 0.693 |
| betA | gene | 0.584 | 0.854 | 0.997 | 0.747 | 0.852 | 0.854 | 0.621 |
| betB | gene | 0.487 | 0.821 | 0.923 | 0.74 | 0.82 | 0.821 | 0.538 |
| gshA | gene | 0.77 | 0.971 | 0.997 | 0.947 | 0.97 | 0.971 | 0.837 |
| gshF | gene | 0.664 | 0.974 | 0.977 | 0.97 | 0.973 | 0.974 | 0.79 |
| gshB | gene | 0.898 | 0.964 | 0.989 | 0.941 | 0.961 | 0.964 | 0.913 |
| gor | gene | 0.955 | 0.984 | 0.993 | 0.975 | 0.98 | 0.984 | 0.961 |
| cspA | gene | 0.778 | 0.884 | 0.992 | 0.797 | 0.881 | 0.884 | 0.813 |
| trkA | gene | 0.804 | 0.873 | 0.817 | 0.937 | 0.87 | 0.873 | 0.819 |
| ktrA | gene | 0.778 | 0.837 | 0.82 | 0.856 | 0.833 | 0.837 | 0.791 |
| trkH | none | 0.101 | 0.11 | 0.06 | 0.696 | 0.11 | 0.11 | 0.103 |
| ktrB | none | 0.077 | 0.171 | 0.112 | 0.359 | 0.173 | 0.171 | 0.086 |
| ktrD | none | 0.144 | 0.202 | 0.126 | 0.512 | 0.204 | 0.202 | 0.165 |

## Findings

1. **Catalogued strains gain most where translated search was weakest.**
   ectB 0.70 -> 0.98, gshF 0.66 -> 0.97, betB 0.49 -> 0.82, betA 0.58 ->
   0.85, mrpE 0.69 -> 0.93. With A, 1 reported unit stays below 0.5
   (betL), against 2 without a catalog (proP, betB).
2. **betL is the one regression** (0.868 -> 0.459). The catalog inherits the
   full-length HMM's decision, and for betL that HMM cannot separate betL
   from betT/caiT (annotate HMM F1 0.55). Translated search gets betL right
   through DIAMOND's betL_decoy references.
3. **Five excluded units become reportable from catalog reads**
   (catalog only, A): murB 0.97, mrpA 0.85, otsA 0.69, mscS 0.64, mrpC 0.57.
   In `unit_counts.tsv` a unit on none already reports its catalog reads
   (status `catalog_only`), so this needs no new mechanism.
4. **Some units stay unreliable even at full length:** opuAC/opuBC/opuCC
   (<= 0.09), the opu ATPases (0.32-0.37), mrpD (0.39) and the trk/ktr
   channel (0.11-0.20). Their paralogs confuse the HMMs at full length too.
   They stay excluded.
5. **The release's small index is as good as the full one here** (0.892
   both). Its only cost is background false positives (270 vs 170 per
   million): reads from genes without an HMM hit are not in the small index,
   so they fall back to translated search. In a real catalog, where almost
   every gene has no hit, that gap will be larger than in this simulation.
6. **The fallback barely matters when strains are catalogued** (catalog only
   0.890 vs hybrid 0.892). Its role is the uncatalogued strains (B).
7. **B settles more than expected** (107,302 test reads, 28,743 of them
   positives, 99.6% correctly called). They are near full-length alignments
   at >= 97% identity to training genes of the same family. CD-HIT at 90%
   limits identity over the whole protein, not over a 42-residue stretch.
   Because the simulation uses one codon per amino acid, nucleotide identity
   equals protein identity here. Real strains differ synonymously and would
   reach 97% far less often, so B's gain over no catalog (0.800 vs 0.770) is
   optimistic. The realistic lower bound is the no-catalog line.

## Two bugs this benchmark found (fixed in `df87119`)

The first run (kept in `results/catalog_benchmark_v1/`) gave A a median of
0.871, with mscL, kdpA and nhaA worse than without a catalog. Both causes
were in osmotool's catalog assignment, not in the approach:

- **Gathering cutoff.** A gene needed one domain at or above GA2 as well as
  its sequence score above GA1. `hmmsearch --cut_ga`, which annotate uses,
  reports a sequence at GA1 even when its hit is split into domains each
  below the cutoff. Example: a kdpA scoring 743 bits against GA 530, as
  domains of 249 + 499 bits. Now the sequence score alone decides, as in
  annotate.
- **Envelope.** Only reads inside the HMM envelope counted. The models are
  built from trimmed alignments, so the envelope misses the gene ends
  (nhaA: 100-846 of a 1359 nt gene). Now a read anywhere on an annotated
  gene counts, as in translated search. On a fused ORF where a second family
  passes too, a read in that family's region goes to it.

## Caveats

- A is an upper bound: the catalog holds the very strain each read came
  from. Real catalogs hold relatives at various distances. Only the two
  ends are scored, because the single-codon simulation cannot represent a
  related strain's nucleotide divergence.
- A real catalog is mostly genes without an HMM hit. Here, 83% of the A
  catalog's genes have a hit, by construction. The background test reads
  cover the unrelated-gene case only through the 4,992 random background
  proteins.
- Catalog calls are full-length HMM decisions, so this benchmark is also a
  read-level view of annotate's errors. Hard-negative calls are scored
  against UniProt-derived labels, which are known to be imperfect for some
  families.

## Update: betL decoy (HMM decoy built from betL's negatives)

Open decision 1 below, tried: `betTcaiT_decoy`, an HMM built from
`betL.negative.train` (every UniProt betT/caiT entry is in it), with betL's GA
calibrated decoy-aware (`families.yaml` `hmm_decoys`, `05b`, `06`). Test
release `releases/v10-betl-decoy`: the same refs, splits and alignments as v10,
rebuilt and benchmarked; a catalogue benchmark rerun.

| betL, f1_all (precision_all / recall) | v10 | with decoy |
|---|---|---|
| profile, no catalog (DIAMOND + cascade) | 0.868 (0.929 / 0.815) | unchanged |
| profile, HMM only | 0.584 (0.458 / 0.805) | 0.740 (0.925 / 0.617) |
| annotate, HMM (protein level, F1) | 0.554 (0.547 / 0.562) | 0.841 (1.000 / 0.726) |
| catalog A (strains in the catalog) | 0.459 (0.465 / 0.454) | **0.857** (1.000 / 0.750) |
| catalog B (not in the catalog) | 0.776 (0.802 / 0.752) | 0.852 (0.932 / 0.785) |

- betL's GA drops from 570.8 to 369.1 bits. The decoy takes 17 of 73 calibration
  positives (betL proteins that score higher against betT/caiT than against
  betL's own model) and 377 of 384 negatives. The calibration F1 of 0.82 is
  optimistic, since the decoy was trained on those negatives; the held-out test
  numbers above are the real ones.
- Nothing else moves: annotate HMM F1 is unchanged for the other 43 families
  (overall micro F1 0.928 -> 0.931), the cascade config and the report levels
  are identical to v10, and every other unit's f1_all is unchanged in both the
  profile and the catalogue benchmark. Median f1_all over the 27 reported units:
  0.770 / 0.892 / 0.800 (no catalog / A / B) in v10 and 0.770 / 0.892 / 0.804
  with the decoy.
- Caveats: betL has only 95 held-out hard-negative sequences, so a precision of
  1.000 rests on a small sample. Recall is capped at about 0.73-0.75 (annotate
  and catalogue A), because the decoy takes some real betL proteins; with HMM
  alone in profile mode recall even falls (0.805 -> 0.617), though precision
  rises more. betL is one of the NaCl model's 46 features: adopting the decoy
  changes betL counts for every annotated genome, so the BacDive annotation and
  the NaCl model would have to be redone.

## Open decisions

1. **betL in the catalog.** Options:
   - (a) annotate the catalog per family with the method that is better at
     full length on calibration proteins (DIAMOND with decoys for betL);
   - (b) treat betL catalog reads as fallback, so they go through DIAMOND;
   - (c) accept it and document it.

   (a) is the general fix.
2. **Report levels with a catalog.** No recalibration seems needed for the
   fallback path. The reads it sees in practice come from uncatalogued
   strains, which is what v10's calibration reads represent, and B shows the
   translated numbers hold. Units on none could be shown as reportable from
   catalog reads (finding 3); `unit_counts.tsv` already gives their catalog
   counts.
3. **Defaults.** `--catalog` is not the default yet. Remaining questions:
   - which catalog(s) a release ships;
   - whether profile uses the release's catalog when it has one;
   - DIAMOND's thresholds for the fallback (osmotool defaults to 0.80/0.80,
     the benchmark runs 0.40/0.50).
