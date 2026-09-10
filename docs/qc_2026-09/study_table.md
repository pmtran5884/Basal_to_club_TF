## Studies contributing to this run

| study (HLCA label) | donors basal | donors club | donors goblet | donors ciliated | secretory cells scored | atlas club ref cells | atlas goblet ref cells | club/goblet AUROC | used for |
|---|---|---|---|---|---|---|---|---|---|
| Barbry_Leroy_2020 | 9 | 10 | 4 | 9 | 9,972 | 7,095 | 1,116 | 0.451 | airway subset; secretory resolver; pseudobulk:basal; pseudobulk:club; pseudobulk:goblet; pseudobulk:ciliated; LOSO eval |
| Banovich_Kropski_2020 | 2 | 17 | 0 | 24 | 3,959 | 813 | 151 | 0.687 | airway subset; secretory resolver; pseudobulk:basal; pseudobulk:club; pseudobulk:ciliated; LOSO eval |
| Seibold_2020 | 14 | 6 | 6 | 13 | 3,189 | 496 | 2,189 | 0.578 | airway subset; secretory resolver; pseudobulk:basal; pseudobulk:club; pseudobulk:goblet; pseudobulk:ciliated; LOSO eval |
| Nawijn_2021 | 10 | 6 | 6 | 10 | 2,838 | 1,782 | 651 | 0.588 | airway subset; secretory resolver; pseudobulk:basal; pseudobulk:club; pseudobulk:goblet; pseudobulk:ciliated; LOSO eval |
| Misharin_2021 | 1 | 2 | 1 | 2 | 1,933 | 23 | 1,112 | 0.825 | airway subset; secretory resolver; pseudobulk:basal; pseudobulk:club; pseudobulk:goblet; pseudobulk:ciliated; LOSO eval |
| Krasnow_2020 | 3 | 3 | 1 | 3 | 1,506 | 51 | 22 | 0.765 | airway subset; secretory resolver; pseudobulk:basal; pseudobulk:club; pseudobulk:goblet; pseudobulk:ciliated; LOSO eval |
| Misharin_Budinger_2018 | 0 | 1 | 0 | 6 | 437 | 4 | 9 | — | airway subset; secretory resolver; pseudobulk:club; pseudobulk:ciliated; LOSO eval: too few ref cells |
| Lafyatis_Rojas_2019 | 0 | 1 | 0 | 5 | 314 | 122 | 21 | 0.585 | airway subset; secretory resolver; pseudobulk:club; pseudobulk:ciliated; LOSO eval |
| Teichmann_Meyer_2019 | 0 | 0 | 0 | 2 | 73 | 4 | 11 | — | airway subset; secretory resolver; pseudobulk:ciliated; LOSO eval: too few ref cells |
| Meyer_2019 | 0 | 0 | 0 | 2 | 64 | 5 | 0 | — | airway subset; secretory resolver; pseudobulk:ciliated; LOSO eval: too few ref cells |

## Datasets specified in `config/datasets.yaml`

| dataset id | role in design | source | accession / dataset id | note |
|---|---|---|---|---|
| deprez_2020 | atlas | geo | GSE143868 | Positional series. airway_position is real biology and must enter the model as a covariate; nasal locations are flagged, never silently pooled with bronchial. |
| goldfarbmuren_2020 | atlas | geo | GSE134174 | — |
| ruizgarcia_2019 | atlas | geo | GSE121600 | Ground truth for both the differentiation trajectory (real ALI time, not pseudotime) and the club/goblet resolver. GSM3439925 (in vivo bronchial biopsy) is the designated resolver holdout. Media is a covariate with a mandatory per-medium sensitivity run: BEGM and PneumaCult give different secretory composition. |
| ravindra_2021 | atlas | cellxgene | GSE166766 | The genuinely bronchial ALI set. Mock/uninfected cells ONLY: infected cells carry an interferon program that will dominate any differential signature. The infection filter is asserted at load time, not left to a downstream subset. |
| hlca_core | reference | cellxgene | 066943a2-fdac-4b29-b348-40cede398e4e | scArches/scANVI label-transfer anchor. Never merged into the atlas as a batch. |
