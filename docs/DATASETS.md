# Candidate datasets — `Basal_to_club_TF`

Retrieved from PubMed→GEO links, a GEO DataSets screen (234 human airway scRNA-seq series), and the CELLxGENE Discover collection index (389 collections → 127 human lung/airway datasets). Every accession, sample count and cell count below came back from those queries; nothing here is quoted from memory. Machine-readable version: `candidate_datasets.csv`.

**How to use this page:** pick a set from Tier A + Tier B, decide on Tier C, and I'll write the loaders. My recommended minimum is marked ✅.

---

## Ground truth — your pick, confirmed as a good one

### `GSE121600` — Ruiz García et al., *Development* 2019 (PMID 31558434) ✅
"Single-cell RNA-seq reveals novel cell differentiation dynamics during human airway epithelium regeneration"

This is a better fit than the title suggests. The series is not just one culture — its stated design is *"human in vitro time course differentiation, human bronchial biopsy, human nasal brushing, human nasal turbinate, mouse in vitro differentiation and pig trachea"*. The 22 samples break down as:

| Sample block | GSMs | Content |
|---|---|---|
| BEGM ALI time course | GSM3439913–3439919 | ALI day 2, 4, 7, 12, 17, 22, 62 |
| PneumaCult ALI time course | GSM3439920–3439923 | ALI day 7, 12, 28, 47 |
| Additional donors | GSM4080085–4080090 | BEGM ALI d14 ×3 (donors D252, D266, D267), BEGM ALI d31 (D395), PneumaCult ALI d26 (D395), d28 (D389) |
| Matched in vivo | GSM3439925–3439927 | bronchial biopsy, nasal brushing, nasal turbinate |
| Non-human | GSM3439924, GSM3439928 | mouse tracheal ALI d3, pig trachea — **excluded** |

Why it works as ground truth for the parts you named:
- **Trajectory ground truth.** A dense ALI series from day 2 (basal-dominated) to day 62 (fully differentiated) is a directly observed basal→secretory/ciliated transition, not a pseudotime inferred from a steady-state snapshot. Real time is the strongest available anchor for "drives".
- **Resolver ground truth.** The secretory compartment is sampled repeatedly across two media, so the club/goblet split can be validated against a population whose composition shifts for a known reason rather than against a single set of author labels.
- **In vitro ↔ in vivo bridge.** Bronchial biopsy, nasal brushing and turbinate in the *same series and pipeline* let you test whether the culture-derived club program matches the in vivo one without a cross-study batch effect confounding the comparison.

Two things to handle explicitly, both of which become config, not caveats:
1. **Filter to human.** The mouse and pig samples are in the same accession.
2. **Media is a real covariate, not noise.** BEGM and PneumaCult produce different secretory compositions. Treat `media` as a batch covariate and run the club prediction separately per medium as a sensitivity check — if the club TF list is medium-dependent, that is a finding you need to know about before it reaches a figure.
3. The cultures are airway/nasal-derived rather than strictly bronchial. Pairing this with a genuinely bronchial ALI set (Ravindra, below) covers that gap.

---

## Tier A — in vivo backbone (steady-state composition, donor diversity)

| ✅ | ID | Accession | Study | Size | Why it's here | Watch out for |
|---|---|---|---|---|---|---|
| ✅ | `deprez_2020` | **GSE143868** | Deprez et al., *AJRCCM* 2020 (PMID 32726565) — "Single cell RNA-seq mapping of nasal and tracheobronchial airways in human healthy volunteers" | 35 samples, 77,969 cells, 10 healthy volunteers, 35 locations nose→12th division | The strongest single in vivo backbone: healthy-only, positionally resolved, and cold-active protease dissociation gave 89.1% epithelial cells — you get secretory and ciliated cells in proportion rather than a basal-biased survivor population | The nose→bronchus gradient is real biology; it must be a covariate, and nasal locations should be flagged rather than silently pooled |
| ✅ | `goldfarbmuren_2020` | **GSE134174** | Goldfarbmuren et al., *Nat Commun* 2020 (PMID 32427931) — smoking effects and lineage reconstruction in human tracheal epithelium | 28 samples | Purpose-built for lineage reconstruction from basal cells, large donor count, and the nonsmoker arm is a clean healthy tracheal set | Smoker/nonsmoker split — nonsmokers primary, smoking as covariate or a sensitivity arm |
|  | `madissoon_2023` | CELLxGENE `0e9d47fb-89b1-42d8-b426-2c7630b5f5fa` | "A spatially resolved atlas of the human lung characterizes a gland-associated immune niche" (doi 10.1038/s41588-022-01243-4), *Airway epithelial cells* subset | 29,505 cells | Curated h5ad with ontology-standardized labels, and it explicitly resolves **submucosal gland** cells — the exact contaminant that corrupts a club call (gland mucous cells are MUC5B-high, serous are LTF/LYZ-high) | Take the airway-epithelium dataset only, not the whole multi-tissue collection |
|  | `travaglini_2020` | CELLxGENE `8c42cfd0-…` (10X), `e04daea4-…` (Smart-seq2) | Travaglini et al., *Nature* 2020 (doi 10.1038/s41586-020-2922-4) | 65,662 + 9,409 cells | The Smart-seq2 arm has deep per-cell coverage, which matters specifically for TFs with low mRNA detection — the population VIPER exists to rescue | Whole-lung; proximal airway is a subset. Two chemistries = two batches |

**Annotation reference (not a batch to merge):**

| ✅ | ID | Accession | Notes |
|---|---|---|---|
| ✅ | `hlca_core` | CELLxGENE `066943a2-fdac-4b29-b348-40cede398e4e` | HLCA **core**, 584,944 cells (Sikkema et al., *Nat Med* 2023, doi 10.1038/s41591-023-02327-2). Harmonized multi-level annotation; the natural scArches/scANVI label-transfer anchor |
|  | `hlca_full` | CELLxGENE `9f222629-9e39-47d0-b83f-e08d610c7479` | HLCA **full**, 2,282,447 cells. Use only if the core under-represents a query dataset |

---

## Tier B — in vitro ALI (where the basal→club transition is actually observable)

| ✅ | ID | Accession | Study | Size | Role |
|---|---|---|---|---|---|
| ✅ | `ruizgarcia_2019` | **GSE121600** | see above | 18 human samples | Ground truth + primary trajectory |
| ✅ | `ravindra_2021` | **GSE166766** / CELLxGENE `030faa69-ff79-4d85-8630-7c874a114c19` | Ravindra et al., *PLoS Biol* 2021 (PMID 33730024) | 77,650 cells (curated h5ad) | The genuinely **bronchial** HBEC ALI set, longitudinal, with a mock arm. Also the bridge to your virus work. **Mock/uninfected cells only** — infected cells carry an interferon program that will dominate any differential signature |
|  | `barbry_ali_conditions` | CELLxGENE `c69fb6cd-fc4d-4216-85cb-8d80e7771786` | "Cell culture differentiation and proliferation conditions influence the in vitro regeneration of the human airway epithelium" (doi 10.1165/rcmb.2023-0356MA) | 10,224 cells | Companion to the ground-truth series; quantifies how culture conditions move the secretory compartment. Directly informs the BEGM-vs-PneumaCult covariate |
|  | `ali_organoid_2024` | CELLxGENE `ef1a4b05-540c-4f73-8ceb-49dfa800645d` | "Transcriptomic Analysis of Air-Liquid Interface Culture in Human Lung Organoids Reveals Regulators of Epithelial Differentiation" (doi 10.3390/cells13231991) | 15,806 cells | Independent differentiation system — tests whether club drivers are culture-system-invariant |
|  | `hbec_ali_flu_2026` | GSE319472 | HBEC ALI ± H3N2 (PA-X) | 12 samples | Extra HBEC ALI donors; mock arm only |
|  | `lung_ali_cov2_2025` | GSE272756 | "SARS-CoV-2 infection of human lung ALI cultures reveals basal cells as relevant targets" | 54 samples | Large ALI sample count; mock arm only; verify proximal vs distal origin per sample |

---

## Tier C — disease/perturbation (control arms in, disease arms as external tests)

| ID | Accession | Study | Size | Role |
|---|---|---|---|---|
| `carraro_2021` | GSE150674 | Carraro et al., *Nat Med* 2021 (PMID 33958799) — CF airways at single-cell resolution | 50 samples | Non-CF controls add donors; the CF arm is a held-out test of whether the club program moves as expected in a disease with a remodeled secretory compartment |
| `cf_healthy_biopsy` | CELLxGENE collection *"Cystic Fibrosis and healthy control biopsy"* | — | 96,479 cells (carina of trachea, respiratory tube) | Large curated proximal-airway biopsy set with healthy controls |
| `copd_ali_2026` | GSE302794 | Healthy / ex-smoker COPD / current-smoker COPD HBEC in ALI | 4 samples | Healthy arm is true HBEC ALI; very small n |

---

## Tier D — reference only, never merged into the atlas

| ID | Accession | Use |
|---|---|---|
| `plasschaert_2018` | GSE102580 | Plasschaert et al., *Nature* 2018 (PMID 30069046) — origin of the human pulmonary ionocyte definition. Grounds the ionocyte marker/TF panel and serves as a rare-cell sanity set. Mixed human/mouse series |
| `montoro_2018` | GSE103354 | Montoro et al., *Nature* 2018 (PMID 30069044) — **mouse trachea**. Source for the hillock/tuft/ionocyte lineage-TF definitions only. Not merged; cross-species marker transfer is imperfect |

---

## Recommended minimum (the ✅ set)

`GSE143868` + `GSE134174` (in vivo, healthy, many donors) + `GSE121600` human samples (ALI time course, ground truth) + `GSE166766` mock arm (bronchial ALI) + HLCA core as label-transfer reference.

That is two in vivo sources across nose→bronchus and trachea, two in vitro sources across two culture systems and two media, and one annotation anchor — enough donor and platform diversity for leave-one-dataset-out to mean something, without so many batches that integration becomes the dominant source of variance.

**Explicitly excluded by default:** distal/alveolar-only datasets, nasal-only datasets (except where they arrive inside a proximo-distal series like GSE143868, and then flagged), whole-lung datasets with negligible proximal airway, and every infected/stimulated arm.

---

## Still needed from you

1. Confirm the ✅ set, or edit it.
2. Any **in-house HBEC data** to add as an anchor batch? (Not required — the design works without it.)
3. Held-out set for validating the club/goblet resolver: default is to hold out the `GSE121600` in vivo bronchial biopsy sample, since it is ground truth from the same lab and pipeline as the cultures. Alternative is `madissoon_2023`, whose value is that it separately labels submucosal gland.
4. Access route: GEO supplementary matrices are direct downloads; CELLxGENE datasets are curated `.h5ad` with standardized ontology labels — for the CELLxGENE-hosted ones I'd default to the curated h5ad rather than re-processing from raw.
