# Bronchial ALI time course (GSE233145) — methods

This arm repeats the DREM time-course analysis of `docs/drem_ali/` on *bronchial*
cultures, so that the nasal result can be separated from the tissue of origin.
Everything downstream of the loader is shared code; only the ingestion module and
the sampled day grid differ.

## Dataset

| | nasal arm | bronchial arm |
|---|---|---|
| accession | GSE121600 | GSE233145 |
| tissue | nasal epithelium, ALI | bronchial epithelium (primary HBEC), ALI |
| platform | 10x 3' | Drop-seq |
| donors | D246, D275 | Donor_1, Donor_2 |
| sampled days | 2, 4, 7, 12, 17, 22 | 0, 3, 5, 7, 14, 21, 28 |
| cells after QC | 24,766 | 77,714 |
| median genes/cell | 3,574 | 418 |

GSE233145 is deposited as one per-sample Drop-seq DGE text matrix per GSM plus a
single cell-level metadata table. The metadata has no barcode column, so it
cannot be joined to the matrices cell-by-cell; the deposit's own annotation is
therefore used only as a *composition-level* control (see below), and all labels
used in the analysis are computed here.

## Ingestion and QC (`src/basal_to_club/data/gse233145.py`)

1. Per-sample DGE matrices are read in gene chunks and concatenated on the
   **union** of gene symbols (`fill_value=0`). Each file lists only the genes
   detected in that sample, so an inner join silently deletes markers — *FOXJ1*
   and *MUC5B* were both lost this way in a first pass, which would have made the
   ciliated control unrecoverable by construction.
2. QC filters match the published thresholds: `MIN_COUNTS`, `MIN_GENES`,
   `MAX_PCT_MITO` as recorded in the module.
3. Day and donor are parsed from the GSM titles into `obs.day` / `obs.donor`.

## Labels

Two label arms are carried through the whole analysis, exactly as in the nasal
arm:

* `class_atlas` — HLCA-reference label transfer (`labels/atlas_transfer.py`),
  held-out accuracy 0.881 on this dataset.
* `class_resolver` — the secretory resolver (`secretory/resolver.py`), which
  scores club / goblet / hybrid programs per cell and refuses to call a class
  whose marker panel is depleted. The submucosal-gland programs are optional on
  this dataset: a basal-cell-derived ALI culture cannot contain SMG cells, and
  the panel guard correctly refuses them.

## Controls (pre-DREM)

* **Composition vs the authors' annotation.** Called fractions per sample against
  the deposit's published composition, for the three classes the deposit
  expresses (basal, ciliated, secretory).
* **Differentiation direction.** Class fractions against interface day, per donor.
* **Marker positivity per called class.** *KRT5*, *SCGB1A1*, *MUC5AC*, *MUC5B*,
  *FOXJ1*.
* **Platform depth.** Per-marker detection rate in this dataset against the nasal
  dataset, to bound what the labels can support at Drop-seq depth.

See `bronchial_labeling_controls.png`.

## DREM configuration

Identical to the nasal arm (`scripts/run_drem_ali.py`; patched DREM in
`external/STEM_DREM/patched`): log2 CPM pseudobulk series per
arm x donor x lineage, regulons built from ChEA3 libraries restricted to the
expressed universe and to expressed TFs, two regulon families
(`binding` = ENCODE/ReMap/Literature ChIP-seq, `coexpression` =
ARCHS4/GTEx/Enrichr co-expression), capped at 1000 targets per TF.

**One deviation: the day grid is thinned to days 0, 3, 7, 14, 28**
(`BTC_DREM_DAYS`). DREM's path search grows steeply with the number of sampled
time points: on the full seven-point grid a single `binding` model took 2.4 h
(8.67e6 ms) versus 7.6 min for the corresponding six-point nasal model, which put
the 32-model grid at >30 h. The thinned five-point grid restores nasal-like cost
(392 s for the same model) and keeps the shape of the course, with the dense early
sampling represented by day 3 and day 7.

## Comparison to the static VIPER arm

`scripts/compare_drem_viper.py --drem-set drem_bronchial` — the definitions
(rank percentile, control TF panels, club-specificity rule, reproducibility
statistics, convergence test) are fixed in that script's docstring and are the
same ones used for the nasal arm, which it reproduces exactly when run with
`--drem-set drem_ali`.
