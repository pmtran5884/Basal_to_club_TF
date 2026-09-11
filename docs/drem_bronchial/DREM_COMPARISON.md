# Bronchial time course vs the static VIPER arm

Independent replication of `docs/drem_ali/DREM_COMPARISON.md` on *bronchial* ALI
cultures (GSE233145, primary HBEC, Drop-seq, donors Donor_1/Donor_2, days
0/3/7/14/28 of the sampled 0-28 course). Methods and the one deviation from the
nasal arm (thinned day grid) are in `METHODS.md`. All 32 models
(2 label arms x 2 donors x 4 lineages x 2 regulon families) exited cleanly.

Reproduce with:

```bash
BTC_DREM_SET=drem_bronchial BTC_DREM_DAYS=0,3,7,14,28 python scripts/prep_drem_bronchial.py
BTC_DREM_SET=drem_bronchial BTC_DREM_DAYS=0,3,7,14,28 python scripts/run_drem_ali.py
python scripts/compare_drem_viper.py --drem-set drem_bronchial --label "bronchial GSE233145"
python scripts/plot_drem_cross_dataset.py
```

## 1. The labeling machinery holds up on bronchial cultures

`bronchial_labeling_controls.png`

* **Composition vs the authors' annotation** (15 samples x 3 classes the deposit
  expresses): basal r = 0.92, ciliated r = 0.96, secretory r = 0.82; mean absolute
  error 0.057-0.066 of a sample. HLCA transfer held-out accuracy 0.881.
* **Direction of differentiation**: basal falls (0.83 -> 0.47 / 0.62 by day 28),
  ciliated rises from ~0 to 0.36 / 0.20, in both donors.
* **Marker positivity by called class**: club calls are 80.9% *SCGB1A1*+ against
  25.9% in basal; ciliated calls are 14.3% *FOXJ1*+ against 0.9% in basal.
* **Goblet calls are not supported at this depth**: 0.4% *MUC5AC*+, 0.6% *MUC5B*+.
  The dataset itself is the limit — *MUC5AC* is detected in 0.24% of all cells here
  versus 7.7% in the nasal 10x set, at 418 versus 3,574 median genes per cell. The
  highest mucin positivity in the object sits in the resolver's `hybrid` class
  (1.9% / 4.4%), i.e. what mucous signal exists is being caught, there is just very
  little of it.

**Read this as: basal, club and ciliated are sound on this dataset; goblet is
weak, and that weakness is platform depth, not the resolver.**

## 2. Control TFs: the ciliated arm replicates, the goblet arm does not

Best rank percentile across the four models per lineage (coexpression family;
`tables/drem_control_recovery.csv`):

| lineage | TF | nasal | bronchial |
|---|---|---|---|
| ciliated | *FOXJ1* | 0.002 | **0.003** |
| ciliated | *RFX2* | 0.007 | **0.005** |
| ciliated | *RFX3* | 0.020 | **0.006** |
| ciliated | *TP73* | 0.046 | **0.025** |
| ciliated | *MYB* | 0.527 | **0.025** |
| goblet | *SPDEF* | 0.041 | 0.426 |
| goblet | *XBP1* | 0.078 | 0.694 |
| goblet | *CREB3L1* | 0.440 | 0.396 |
| club | *ELF3* | 0.031 | 0.070 |
| club | *KLF5* | 0.013 | 0.465 |
| club | *CEBPB* | 0.252 | 0.428 |
| club | *SOX2* | 0.574 | 0.517 |

*FOXJ1* is rank 1-3 of ~800 TFs in **every** ciliated model in both tissues, and
the rest of the ciliated panel is tighter on bronchial than on nasal. The goblet
positive control fails here (*SPDEF* at the 43rd percentile), exactly as the
marker-positivity control predicts — with 0.4% *MUC5AC*+ cells there is no mucous
trajectory for DREM to fit. The club priors are mid-pack in both tissues, which is
the expected behaviour of a weak prior rather than a failure.

## 3. Reproducibility improves, but donor still dominates

Mean Spearman of club TF rankings (`tables/drem_club_reproducibility.csv`):

| comparison | nasal | bronchial |
|---|---|---|
| same donor, both label arms | 0.63 (11/20 shared) | **0.83 (17/20)** |
| different donors | 0.23 (6.5/20) | **0.49 (7/20)** |

Swapping the entire labeling arm perturbs the club ranking far less than swapping
donor, in both tissues. The bronchial arm is the more stable of the two — deeper
series (5 time points spanning day 0 to 28, ~78k cells) rather than better labels.

## 4. The time-course evidence still does not corroborate the VIPER club set

* Of the 11 robust club-specific TFs from the static VIPER arm, 6 are scored in the
  bronchial club models. Their median club rank percentile is 0.492 against 0.501
  for the other 3,076 TF-model observations — one-sided Mann-Whitney p = 0.28
  (nasal: 8 of 11, 0.480 vs 0.500, p = 0.30).
* The strict club-specificity rule (top decile in the club model, outside the top
  quartile of both the goblet and the ciliated model, in >=3 of 4 arm x donor
  combinations) returns **zero** TFs on bronchial, as it did on nasal.
* Across tissues the club rankings themselves agree only weakly: Spearman 0.37 over
  570 shared TFs, 2 of 20 top-20 TFs shared
  (`tables/cross_dataset_club_agreement.json`).

`drem_controls_nasal_vs_bronchial.png` shows all four of these comparisons.

**The only TFs in the top decile of the club models in both tissues are *FOXM1*,
*MYBL2*, *TFDP1*, *PAX9* and *ZBTB7C*.** Three of those five (*FOXM1*, *MYBL2*,
*TFDP1*) are cell-cycle/proliferation regulators, i.e. the reproducible part of the
basal->club time-course signal is proliferation, not club identity. *PAX9* and
*ZBTB7C* are the only two candidates that are both reproducible across tissues and
not obviously cell-cycle.

A further caveat specific to this dataset: the bronchial club models are topped by
*FOXJ1*, *RFX2*, *RFX3* and *CCDC17* — ciliated regulators. At 418 genes per cell
the club calls carry ciliated signal (early/transitional ciliated cells are hard to
separate from club cells at this depth), so the bronchial club branch is partly a
ciliated branch. This is a second reason to treat the club-arm output of the
time-course method as unresolved rather than negative.

## 5. Conclusion

Two independent airway ALI time courses, different tissue and different platform,
agree on the same three points:

1. The pipeline's labeling and control machinery works — the ciliated positive
   control is recovered at rank 1-3 of ~800 TFs everywhere it can be tested, and
   composition transfer reproduces the authors' own annotation.
2. Whether the goblet positive control is recovered is set by sequencing depth
   (nasal 10x yes, bronchial Drop-seq no), not by the labeling strategy.
3. Neither time course supports the static VIPER club-specific TF set. VIPER's club
   TFs are median-ranked in the DREM club models of both tissues, and the DREM club
   rankings do not agree with each other across tissues beyond a proliferation
   signature.

The practical consequence for the repository: the VIPER club ranking should be
reported as a hypothesis list that the orthogonal time-course method does not
confirm, and the club arm needs an experiment that separates club from
transitional-ciliated identity at higher depth (10x or deeper Drop-seq, ideally
with a club-lineage reporter) before any candidate is promoted.
