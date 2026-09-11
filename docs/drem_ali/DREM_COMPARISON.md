# DREM on ALI time courses vs VIPER on the static atlas

Second, independent line of evidence for the basal->club regulators. VIPER scores a
TF on the *state* difference between basal and club cells in a static atlas; DREM
(Dynamic Regulatory Events Miner) fits a branching hidden Markov model to the
*time course* of an ALI culture and asks which TFs' targets are over-represented on
each ascending branch. The two use the same regulon families, so a TF that only
appears in one of them is a method-dependent call, not a consensus driver.

## Design

* Data: GSE121600 nasal ALI differentiation time course, donors D246 and D275,
  16 day-series, 24,766 cells, 10,110 expressed genes.
* Labels: cells were assigned to classes by transfer from the HLCA core
  (regularised multinomial logistic classifier, held-out accuracy in
  `results/drem_ali/prep_summary.json`), then pseudobulked per day to give DREM
  its time series. Two label arms, exactly the definitions used by the VIPER
  class-definition arms (`docs/atlas_label_arms/ARM_COMPARISON.md`):
  **atlas label** (club/goblet from the atlas finest labels) and
  **resolver** (the transferred secretory compartment re-split by
  `basal_to_club.secretory.resolve_secretory`).
* Lineages: `all`, `ciliated`, `club`, `goblet` (basal cells plus that lineage only).
* Regulon families: ChEA3 `binding` (ENCODE/ReMap/Literature ChIP-seq) and
  `coexpression` (ARCHS4/GTEx), restricted to TFs expressed in the culture --
  210 and 613 TFs respectively. An unexpressed TF cannot be driving the culture
  and every extra TF is one more covariate in DREM's per-split regression.
* 2 arms x 2 donors x 4 lineages x 2 families = 32 models, all fitted
  (`tables/drem_runs.csv`); 30 produced path (edge) tables, 1 fell back to split
  tables, 1 is a refit.

## Result 1 -- the ciliated control is recovered crisply; goblet partially

Rank percentile within each lineage model (lower is stronger), coexpression
regulons, best over the four arm x donor models (`tables/drem_control_recovery.csv`):

| lineage | TF | best percentile |
|---|---|---|
| ciliated | FOXJ1 | 0.002 (rank 1 in all four models) |
| ciliated | RFX2 | 0.007 |
| ciliated | RFX3 | 0.020 |
| ciliated | TP73 | 0.046 |
| goblet | SPDEF | 0.041 (D246; 0.28-0.57 in D275) |
| goblet | XBP1 | 0.078 |
| club | KLF5 | 0.013 |
| club | ELF3 | 0.031 |

FOXJ1 is the top-ranked regulator of the ascending ciliated branch in every
model, and the rest of the ciliogenesis module follows it inside the top few
percent -- the positive control the pipeline was asked to pass, passed on
time-course evidence rather than state evidence. SPDEF behaves as designed in
D246, the donor with a substantial goblet population, and washes out in D275,
where the atlas goblet calls are poorly supported by mucin expression
(`docs/drem_ali/tables/` and the marker-positivity table in
`results/drem_ali/prep_summary.json`).

Caveat: FOXJ1, RFX3, TP73 and SPDEF are absent from the ChEA3 binding libraries,
so the binding family cannot be scored on these controls at all. Only the
coexpression family carries the controls; binding results are reported for
completeness.

## Result 2 -- the club branch is donor-driven, not definition-driven

Agreement of the club-branch TF ranking (`tables/drem_club_reproducibility.csv`):

| comparison | mean Spearman | mean top-20 overlap |
|---|---|---|
| same donor, different label arm | 0.63 | 11/20 |
| different donor | 0.23 | 6/20 |

This is the opposite pattern to VIPER, where the label definition changed the
club-specific set almost entirely while leaving the global ranking intact. Here
the definition is nearly irrelevant within a donor (D246: 19/20 top-20 shared
across arms) and the donor dominates. No TF is club-specific -- better than the
10th percentile on the club branch and worse than the 25th on both the goblet and
ciliated branches -- in 3 of 4 arm x donor models; 13 TFs manage it in one.

The top club-branch hits are dominated by squamous/basal-differentiation factors
(GRHL1/2/3, ZNF750, IRF6, OVOL1/2, KLF5, EHF, BARX2, PAX9). In these nasal ALI
cultures that program is real and large -- it is also what the atlas transfer
assigns to the hillock-like/suprabasal compartment -- and it is shared between
the club and goblet branches, which is why nothing separates them.

## Result 3 -- the two methods do not converge on club regulators

The 11 definition-robust club-specific TFs from VIPER (8 of which are in the DREM
coexpression TF universe) are not enriched among strong club-branch regulators:
median percentile 0.48 vs 0.50 background, Mann-Whitney p = 0.30.

Read together with Result 1, the honest reading is that the ciliated program is
recoverable by either method, and the basal->club call is not yet
method-independent. The consensus list should be restricted to TFs with
state evidence (VIPER) *and* time-course evidence, and on present data that
intersection is empty at any defensible threshold.

## Patches to DREM

DREM's batch writer crashes twice on these models; both are null dereferences on
state the interactive viewer fills in and batch mode does not. The single-class
source patch is in `external_patches/drem_batch_orderA_nullguard.patch` and is
compiled ahead of the shipped jar on the classpath (`basal_to_club.drem.run`,
`classes_dir=`):

1. `theTreeNode.orderA` (child ordering) is null in batch mode; a local ordering
   by node mean is substituted. Without this, no path table is ever written.
2. `theTreeNode.dpvals` (per-node regulator p-values) is null at nodes where no
   regression was fit; the split-table writers dereference it. The table is now
   skipped for such a node and the recursion continues.

`basal_to_club.drem.parse` prefers path tables and falls back to split tables, so
a model that crashes partway still yields a ranking from whatever DREM wrote.

## Reproducing

    python scripts/prep_drem_ali.py        # transfer labels, pseudobulk per day
    python scripts/run_drem_ali.py         # 32 models, resumes on existing output
    python scripts/collect_drem_ali.py     # rebuild rankings without refitting
    python scripts/refit_drem_run.py atlas_label D275 goblet binding
