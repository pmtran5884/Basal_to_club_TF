# Run report — HLCA core, prior-only arm

Executed 2026-09-09 on a 10-core laptop. This is a **scoped run**, not the full
workflow: it enters the DAG below integration and omits the de-novo (ARACNe)
network arm. Every number below comes from the files in this directory.

## What ran

| Stage | Status | Detail |
|---|---|---|
| Ingest | done | HLCA core, `data/raw/hlca_core/hlca_core.h5ad`, 584,944 cells |
| Airway subset | done | 137,413 cells, 27,402 genes, 10 studies |
| Integration | **skipped** | HLCA is already scVI-integrated; not re-run |
| Secretory resolver | done | 24,285 secretory-compartment cells re-partitioned |
| Pseudobulk + limma-voom | done | donor-blocked, 3 contrasts |
| Prior regulons (ChEA3) | done | 6 libraries |
| VIPER + metaVIPER | done | 1,534 TFs, 3 contrasts |
| De-novo ARACNe arm | **not run** | no ARACNe jar on this machine |
| Benchmark gate | done | **FAILED** — see below |
| Club prediction | done, **provisional** | reported because the gate's failure modes are diagnosable, not because it passed |

## Secretory resolver

The atlas collapses club and goblet into "Secretory" at level 3 and splits the
finest level by anatomical site (`Club (nasal)`, `Goblet (bronchial)`, …) rather
than by secretory program — exactly the situation the resolver exists for.

Resulting classes: basal 81,062 · ciliated 31,344 · **club 9,880** ·
hybrid 4,861 · goblet 1,714 · secretory_unassigned 7,372.

Notable reassignments (`resolver_vs_atlas_label.csv`): 3,246 cells the atlas
calls `pre-TB secretory` are called club; 1,132 `Goblet (bronchial)` cells are
called hybrid rather than either endpoint; 268 `Goblet (subsegmental)` cells are
flagged as ciliated contaminants and removed from the goblet contrast.

## Benchmark gate: FAILED

Primary controls, expression-derived ChEA3 networks (of 1,534 TFs):

| Contrast | TF | Rank | NES | FDR | Verdict |
|---|---|---|---|---|---|
| basal → ciliated | FOXJ1 | **1** | +30.8 | 6.7e-92 | pass |
| basal → goblet | SPDEF | **43** | +3.93 | 0.024 | fail (gate requires ≤ 25) |

Secondaries in the top 10%: goblet 3/5 (XBP1 #4, ELF3 #5, FOXA3 #17, CREB3L1 #32;
EHF #695 fails) · ciliated 3/6 (RFX3 #3, RFX2 #5, TP73 #34; MYB #186, E2F4 #1279
fail, GMNC absent from ChEA3).

Three distinct failure modes, none of them "the biology is absent":

1. **SPDEF at rank 43** clears the 5% percentile threshold (2.8%) but not the
   absolute-rank threshold of 25. Combined across all six libraries it is rank
   58/1,534.
2. **Cross-family confirmation is unmeetable here.** SPDEF and FOXJ1 have no
   regulon in *any* of the three ChIP-seq-derived ChEA3 libraries (ENCODE,
   Literature, ReMap), which cover 336 TFs in total. The frozen requirement
   `require_both_network_families` cannot be satisfied without the de-novo arm.
   This is the concrete argument for the ARACNe arm the design specifies.
3. **Ciliated secondaries fall below 0.5** (0.42), driven by E2F4 scoring
   negative and GMNC being absent from the TF universe.

The discriminative test **passes**: SPDEF ranks 58 in goblet vs 102 in club, and
FOXJ1 ranks 1 in ciliated vs 212 in club. The signal is lineage-specific, not
generic "no longer basal".

`basal_identity` (TP63) and `basal_to_deuterosomal` (MCIDAS) were **not
evaluated** — those contrasts were not computed in this scoped run. MCIDAS is
absent from all six ChEA3 libraries in any case.

## Club prediction (provisional)

`club_tf_ranking.csv`, 1,533 TFs. No TF reaches `high_confidence`, because that
tier requires ≥2 signature constructions and this run used one (endpoint
contrast only): 782 `candidate`, 751 `weak`.

Top 20 by metaVIPER NES: CREB3L4, MITF, XBP1, FOXA1, FOXA2, CREB3L2, MECOM,
ELF3, FOXP4, ZNF189, STAT1, ZNF587, NCOA2, STAT3, SMAD4, DMRTA1, MESP1, NFE2,
ZNF69, HNF4G.

Four of the prior-plausible club TFs land in the top 50 without ever being used
in the ranking: XBP1 (#3), FOXA2 (#5), ELF3 (#8), CREB3L1 (#50). CREB3L4 and
CREB3L2 at #1 and #6 are the same bZIP secretory-program family as CREB3L1.
Counter-evidence worth stating: SOX2 (#1,066), HES1 (#1,298) and CEBPA (#1,205)
score strongly *negative*, and EHF — a plausible secretory factor — is #1,031.

## How to reproduce

```bash
conda env create -f envs/btc-run.yaml   # or see envs/
python scripts/run_hlca_prior_only.py   # writes results/hlca_prior_only/
python scripts/score_hlca_run.py        # gate + club ranking
```

## Bugs this run exposed (all fixed in `src/`, all covered by tests)

1. `limma_voom_signature` built its contrast from alphabetical design-matrix
   column order while accepting `case`/`control` arguments it never used. Since
   "basal" sorts first against club, goblet and ciliated alike, **every**
   signature would have been sign-flipped and every control TF would have ranked
   at the bottom.
2. `run_pyviper` passed the merged signature matrix to aREA. Contrasts are
   filtered independently, so that matrix always contains NaNs, and aREA
   propagates a NaN across the whole sample — the result was an all-NaN activity
   matrix returned **without an error**.
3. `make_pseudobulk` passed an array to `sc.get.aggregate(by=...)`, which expects
   an obs column name.
4. Both R bridges used `pandas2ri.activate()`, removed in rpy2 3.6.
