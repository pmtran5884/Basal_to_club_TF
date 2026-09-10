# `Basal_to_club_TF` — design outline

**Goal.** Aggregate public human bronchial/airway epithelial scRNA-seq into one harmonized reference, build both a de novo (ARACNe) and a prior-knowledge (Ma'ayan-lab / ChEA3-derived) regulon set, and use VIPER to rank transcription factors driving the **basal → club** transition. Correctness is established *before* the club question is asked, by requiring the pipeline to recover known lineage TFs (SPDEF for basal→goblet, FOXJ1 for basal→ciliated) at high rank on the same machinery with the same settings.

This document is the specification for the repository. No code is written yet.

---

## 0. Design principles

1. **Positive controls gate the analysis.** The club prediction is only reported from a configuration that has already passed the SPDEF / FOXJ1 benchmark. The benchmark is a *frozen* module — its thresholds are declared in `config/benchmark.yaml` before any tuning, and the tuning sweep is over pipeline hyperparameters, never over the control TF list.
2. **The club/goblet distinction is a first-class problem, not a labeling detail.** Public airway datasets vary wildly in how they split the secretory compartment; many collapse everything into "secretory" or mislabel club as goblet. A dedicated re-annotation module with an explicit goblet-exclusion rule is required.
3. **Two networks, cross-checked.** De novo ARACNe regulons are airway-specific but composition-sensitive; ChEA3/CollecTRI priors are stable but not airway-specific. Agreement between them is the strongest evidence; disagreement is reported, not hidden.
4. **Every ranked TF carries its evidence.** NES, FDR, leading-edge targets, regulon size, TF expression fraction in each compartment, and which network produced it.
5. **Reproducible from a cold clone.** Snakemake DAG + pinned conda/Docker environments + a dataset manifest of accessions; no manual steps.

---

## 1. Repository layout

```
Basal_to_club_TF/
├── README.md
├── LICENSE                       # MIT (code); data terms per source in docs/DATA.md
├── environment.yml               # python: scanpy, anndata, scvi-tools, pyviper, decoupler, harmonypy, cellrank
├── envs/
│   ├── r-viper.yml               # R: viper, aracne.networks, Seurat (for .rds ingestion)
│   ├── aracne.yml                # ARACNe-AP (java) / ARACNe3
│   └── docker/Dockerfile
├── config/
│   ├── datasets.yaml             # the dataset registry (accessions, format, loader, priors on labels)
│   ├── params.yaml               # QC, integration, metacell, ARACNe, VIPER hyperparameters
│   ├── markers.yaml              # cell-type marker panels incl. club/goblet discriminators
│   ├── benchmark.yaml            # positive-control TFs + acceptance thresholds (frozen)
│   └── regulons.yaml             # network sources, sign-assignment rules, filters
├── workflow/
│   ├── Snakefile
│   └── rules/{ingest,qc,integrate,annotate,secretory,signature,network,viper,benchmark,report}.smk
├── src/airway_viper/
│   ├── ingest/                   # per-dataset loaders → common AnnData schema
│   ├── qc/                       # per-sample QC, doublets, ambient RNA
│   ├── integrate/                # scVI/scANVI + Harmony fallback; HLCA label transfer
│   ├── annotate/                 # hierarchical marker scoring, reference mapping, arbitration
│   ├── secretory/                # club vs goblet vs hybrid resolver (the crux module)
│   ├── trajectory/               # basal→X lineage assignment, pseudotime, metacells
│   ├── signature/                # GES construction (pseudobulk + metacell contrasts)
│   ├── network/                  # ARACNe runner; ChEA3/CollecTRI → regulon converter; sign assignment
│   ├── viper/                    # msVIPER / metaVIPER / NaRnEA wrappers
│   ├── benchmark/                # control-TF rank scoring, config sweep, pass/fail gate
│   ├── report/                   # tables, figures, HTML report
│   └── utils/                    # gene symbol harmonization, TF universe, IO
├── notebooks/                    # exploratory only; not part of the DAG
├── tests/                        # unit tests + a tiny synthetic fixture for CI
├── docs/
│   ├── DATA.md  METHODS.md  DECISIONS.md  RESULTS.md
└── results/                      # gitignored; DVC/zenodo pointer for released outputs
```

---

## 2. Module 1 — Dataset aggregation

**Registry.** `config/datasets.yaml` is the single source of truth. Each entry:

```yaml
- id: <short_key>
  accession: <GEO/EGA/CELLxGENE/HCA id>
  url: <download endpoint>
  format: h5ad | mtx10x | rds | loom
  loader: src/airway_viper/ingest/<id>.py
  tissue: proximal_airway | trachea | bronchus | ALI_culture | nasal
  system: in_vivo_brush | in_vivo_biopsy | ALI | organoid
  chemistry: 10x_v2 | 10x_v3 | dropseq | smartseq2
  donor_key: <obs column>
  celltype_key: <obs column with author labels>
  disease_filter: [healthy]         # restrict; smoking/asthma/COPD as covariate or exclude
  label_granularity: coarse | fine  # does it separate club from goblet?
  trust_labels: true | false        # if false, labels are advisory only
```

**Selection criteria** (to apply when you name the studies):
- Human, proximal airway epithelium (tracheal/bronchial brushings, biopsies, or ALI cultures), untreated/healthy baseline.
- Contains basal **and** a secretory compartment **and** ciliated cells — required for the positive controls.
- Raw counts available (VIPER needs a real expression signature, not a scaled/regressed matrix).
- Cell count ≥ ~2k epithelial cells per study, ≥2 donors.

**Composition target.** Aim for ≥5 studies, ≥15 donors, ≥60k epithelial cells after QC, with both in vivo and ALI represented — ALI cultures are where basal→club differentiation is most clearly captured in time, in vivo brushings give the steady-state distribution. Record the in vivo/ALI split as a covariate; run a leave-one-system-out sensitivity check later.

**The concrete shortlist lives in `docs/DATASETS.md`** — accessions, sample counts, per-study role and per-study caveats, retrieved from GEO and CELLxGENE rather than recalled. Summary of the recommended set: `GSE143868` (Deprez, in vivo nose→bronchus, 10 healthy volunteers) and `GSE134174` (Goldfarbmuren, tracheal, lineage-reconstruction design) as the in vivo backbone; **`GSE121600` (Ruiz García) human samples as the ALI time-course ground truth**; `GSE166766` mock arm (Ravindra) as the genuinely bronchial ALI set; HLCA core as the label-transfer reference. Nasal-only, distal/alveolar-only, and all infected or cytokine-stimulated arms are excluded by default.

**Ground-truth role of `GSE121600`.** Its human block is a dense ALI regeneration series — BEGM days 2/4/7/12/14/17/22/31/62 and PneumaCult days 7/12/26/28/47 across several donors — plus matched bronchial biopsy, nasal brushing and turbinate in the same accession. Real differentiation time is a stronger anchor for "drives" than pseudotime over a steady-state snapshot, and the in vivo samples let the culture-derived club program be checked against tissue without a cross-study batch effect in the way. Two handling rules: filter out the mouse and pig samples that share the accession, and carry `media` (BEGM vs PneumaCult) as a covariate with a per-medium sensitivity run — the two media give different secretory compositions, so a medium-dependent club list is something to discover before it reaches a figure, not after.

**Harmonization.** Gene symbols → a single reference annotation (Ensembl release pinned in `params.yaml`, HGNC alias resolution, drop non-1:1 mappings); counts as integers in `.layers['counts']`; a controlled vocabulary in `.obs` (`dataset_id`, `donor_id`, `sample_id`, `system`, `chemistry`, `author_label`, `disease`). The output of this module is one AnnData with no imputation and no scaling applied.

---

## 3. Module 2 — QC and integration

- Per-sample QC with adaptive (MAD-based) thresholds rather than global cutoffs: counts, genes, % mitochondrial, % ribosomal. Ambient RNA correction (SoupX/CellBender) matters here specifically because `SCGB1A1` and `MUC5AC` are extremely highly expressed secreted transcripts and bleed into basal and ciliated cells — uncorrected ambient signal will corrupt the club/goblet call. Doublet detection per sample (scDblFinder/Solo); basal–secretory doublets otherwise masquerade as transitional cells.
- Non-epithelial cells (immune, endothelial, stromal) removed after a coarse first-pass annotation; PTPRC/EPCAM/COL1A1 gates.
- Integration: **scVI** latent space with `batch_key = dataset_id × donor_id` (or `sample_id`), covariates for chemistry and system. Harmony on PCA as a cheap fallback/sanity comparison. Integration is used for *neighbors/clustering/label transfer only* — VIPER signatures are always computed from the raw counts within a donor, never from corrected values.
- Integration QC: kBET/iLISI batch mixing vs. bio-conservation (cLISI on high-confidence marker-defined cells); report both so over-correction that erases the club/goblet boundary is visible.

---

## 4. Module 3 — Annotation and the club/goblet resolver

This is the module most likely to decide the result, so it gets three independent opinions and an explicit arbitration step.

**Cell-type hierarchy and its lineage-defining TFs.** This table is the content of `config/markers.yaml` and is also the source of the benchmark TF list in §8 — one table, two uses, so a marker panel and its control TF can never drift apart.

| Compartment | Identity markers (for calling the cell) | Lineage-defining TFs (for benchmarking VIPER) |
|---|---|---|
| **Basal** | KRT5, TP63, KRT14, NGFR, DLK2, MIR205HG | **TP63**, KLF5, SOX2, PAX9, TFAP2A |
| **Suprabasal / parabasal** | SERPINB3, NOTCH3, KRT13-low, SCGB1A1-low | KLF4, KLF5, GRHL3, TFAP2A |
| **Hillock / squamous-like** | KRT13, KRT4, SPRR3, S100A2, ECM1 | KLF4, GRHL1, GRHL3, TP63, SOX15 |
| **Club (secretory)** | SCGB1A1, SCGB3A1, SCGB3A2, BPIFA1, CYP2F1, WFDC2, KLK11, TSPAN8 | *(the question — not a control)* prior-plausible set for interpretation only: NKX2-1, FOXA2, SOX2, CEBPA/B/D, ELF3, EHF, XBP1, CREB3L1, HES1, HEY1 |
| **Goblet** | MUC5AC, TFF1, TFF3, S100P, CEACAM5, AGR2 | **SPDEF**, FOXA3, CREB3L1, XBP1, ELF3, EHF |
| **Deuterosomal** | DEUP1, CDC20B, CCNO, HES6, PLK4 | **MCIDAS**, FOXN4, MYB, E2F4, E2F5, GMNC |
| **Ciliated** | FOXJ1, TPPP3, PIFO, CAPS, DNAI1, SNTN, RSPH1 | **FOXJ1**, RFX2, RFX3, TP73, MYB, E2F4, GMNC, MCIDAS |
| **Ionocyte** (rare) | FOXI1, CFTR, ASCL3, ATP6V1G3, ATP6V0D2, TMEM61, CLCNKB | **FOXI1**, ASCL3 |
| **Tuft / brush** (rare) | POU2F3, AVIL, TRPM5, LRMP, SH2D6, BMX | **POU2F3**, ASCL2, SPIB, GFI1B |
| **PNEC** (rare) | CHGA, CALCA, ASCL1, GRP, PCSK1N, SCG2 | **ASCL1**, NEUROD1, INSM1, FOXA2 |
| **Submucosal gland** (exclude, do not analyze) | serous: LTF, LYZ, PRR4, ZG16B; mucous: MUC5B-high, BPIFB2; myoepithelial: TP63⁺ACTA2⁺KRT14⁺ | SOX9, SOX10, MYB — used only to *identify and remove* these cells |

The submucosal-gland row exists because gland mucous cells are MUC5B-high and gland serous cells are LTF/LYZ-high; both will otherwise be swept into the secretory compartment and blur the club definition, particularly in biopsy-derived datasets. Removing them is a prerequisite for the resolver, not an optional refinement.

**Three opinions:**
1. *Author labels* from `celltype_key` — advisory, weighted by `trust_labels`.
2. *Reference mapping* — scANVI/scArches projection onto the HLCA (or whichever reference you designate), giving a transferred label plus an uncertainty score.
3. *Marker-panel scoring* — `scanpy.tl.score_genes` / decoupler over-representation on `config/markers.yaml`, computed per cell and per cluster on Leiden clustering of the integrated space.

**Arbitration.** Per cluster, take the modal call across the three; flag clusters where they disagree or where transfer uncertainty is high, and resolve those by manual inspection of the marker score matrix recorded in `docs/DECISIONS.md`. Cells in low-confidence clusters are excluded from network inference but retained for reporting.

**Club vs goblet vs hybrid — the explicit rule set.** Author labels frequently collapse these, and some datasets label the whole secretory compartment "club" or "secretory" (and, as you noted, sometimes only "ciliated" is broken out). The resolver operates only on cells that passed into the secretory compartment:

- **Club-defining, positive:** SCGB1A1 (high), SCGB3A1, SCGB3A2, BPIFA1/SPLUNC1, CYP2F1, KLK11, WFDC2, MGP (context-dependent), TSPAN8.
- **Goblet-defining, exclusionary:** MUC5AC (high), TFF1, TFF3, S100P, CEACAM5, LTF (submucosal contamination check), AGR2 (high). A cell scoring above threshold on the goblet module is **not** called club, regardless of SCGB1A1.
- **Ambiguous by design:** MUC5B is expressed by both club and goblet cells in human proximal airway and by submucosal glands — it is **excluded** from both discriminant panels and used only as a descriptive covariate. SCGB1A1⁺MUC5AC⁺ cells are assigned to an explicit **club–goblet hybrid / transitional** class and excluded from the club definition used for the driver contrast.
- **Ciliated contamination:** cells scoring on FOXJ1/ciliary modules are removed from the secretory pool even if SCGB1A1⁺ (secretory→ciliated transitional cells and ambient bleed both produce this).

Thresholding uses a two-component Gaussian mixture per marker module fitted within dataset (so scores are comparable across chemistries), with the split point recorded. Output: `club_core`, `goblet_core`, `hybrid`, `secretory_unassigned`, each with a continuous score, plus a per-dataset composition table showing how many author-"goblet" cells were reassigned to club and vice versa — that table goes in the paper/report as a transparency artifact.

**Validation of the resolver itself.** Ground truth is the `GSE121600` human block (default held-out: its in vivo bronchial biopsy sample, with `madissoon_2023` as the alternative because it separately labels submucosal gland). The resolver must reproduce the held-out labels at ≥0.85 balanced accuracy without seeing them. A second, independent check comes free from the ALI design: club abundance must rise monotonically over ALI day in the time course, and the BEGM/PneumaCult comparison must move the goblet fraction in the direction the culture literature reports — a resolver that gets the labels right but the time course wrong is mis-splitting a continuum. Also confirm the `club_core` population is depleted of MUC5AC and enriched for SCGB3A2 relative to `goblet_core` at the pseudobulk level, and that its abundance tracks expectation across in vivo vs ALI.

---

## 5. Module 4 — Defining the transition and building the signature

VIPER needs a **gene expression signature (GES)** — a differential contrast. Three contrasts are computed, and the club prediction is reported on all three so that the answer isn't an artifact of one framing:

1. **Endpoint contrast (primary).** `club_core` vs `basal`, computed as a paired, donor-level pseudobulk difference (limma-voom / DESeq2 moderated t-statistic, donor as blocking factor, dataset as covariate) then converted to a z-scored GES. Paired-by-donor is important: it removes donor and batch effects from the signature without touching the counts.
2. **Trajectory contrast.** Restrict to the basal → suprabasal → club lineage (CellRank2 macrostate/fate-probability or diffusion pseudotime rooted in basal, with the ciliated and goblet branches excluded); GES = correlation of each gene with club fate probability. This captures the *transition* rather than the endpoints and is the more faithful reading of "drives".
3. **Metacell contrast.** SEACells/metacell aggregation within donor to combat dropout, then the same two contrasts at metacell resolution — this is also what feeds ARACNe.

Parallel GESs are built identically for `goblet_core` vs basal and `ciliated` vs basal — these are the benchmark inputs, and building them through the exact same code path is what makes the benchmark meaningful.

---

## 6. Module 5 — Regulon construction (both networks)

### 5a. De novo: ARACNe
- Input: metacell-aggregated, CPM-normalized expression of the airway atlas, restricted to the epithelial compartment; genes filtered to those detected in ≥5–10% of metacells.
- TF universe: human TFs from the Lambert/Vaquerizas curated TF census (~1,600–1,900), intersected with detected genes. Optionally extend to co-factors and chromatin regulators as a separate "signaling/co-TF" run.
- ARACNe-AP (or ARACNe3) with bootstrapping (100 bootstraps, DPI tolerance 0, p ≤ 1e-8), run **per dataset** as well as on the pooled atlas. Per-dataset networks are the substrate for **metaVIPER**, which is exactly the right tool when the correct context-specific network is uncertain.
- Mode of regulation and likelihood from the Spearman correlation between TF and target across metacells; regulons pruned to a size band (e.g. 25–500 targets) to avoid the size bias in NES.

### 5b. Prior knowledge: Ma'ayan-lab and curated resources

ChEA3 ships its underlying gene-set libraries as separate downloadable files — ChIP-seq-derived (ENCODE, ReMap, curated literature ChEA), coexpression-derived (ARCHS4, GTEx), and Enrichr query co-occurrence. The converter targets those downloadable library files rather than the web API, so the build is version-pinned and reproducible offline. Confirm the current library set and release version at build time from the ChEA3 site (`maayanlab.cloud/chea3`) and record it in `results/provenance/`.

- **Converter:** ChEA3 gene-set libraries are unsigned TF→target sets. To become VIPER regulons they need a mode of regulation. Assign sign per edge from the TF–target Spearman correlation *in the airway atlas metacells*, with edges whose |ρ| falls below a floor given likelihood weight 0 (present but uninformative) rather than being dropped. Likelihood = a function of cross-library support (how many of the ChEA3 libraries contain the edge) × |ρ|. This makes the prior network airway-*aware* without making it airway-*derived*.
- **Integrated composite:** build one regulon set per ChEA3 library, plus a composite rank across libraries in the spirit of ChEA3's own orthogonal-integration design (its central claim is that integrating across the distinct library types recovers the correct upstream TF better than any single library does). Practically: keep the per-library networks as **separate metaVIPER inputs** rather than merging them into one undifferentiated edge set.
- **ChEA-KG** (`chea-kg.maayanlab.cloud`), which exposes signed, directed TF→TF edges as a knowledge graph — used not as a primary regulon source but for the **TF–TF layer**: once a club-driver shortlist exists, ChEA-KG gives the signed regulatory relationships among the hits, which is how you tell a master regulator from a downstream effector in the same module.
- **CollecTRI / DoRothEA** (via `decoupler`) as an independent signed, confidence-graded prior — CollecTRI edges already carry sign, so they serve as a check on the sign-assignment heuristic above.
- All prior networks are version-pinned (download date + checksum recorded in `results/provenance/`).

### 5c. Network QC
Report per network: TF count, median regulon size, overlap with the other networks (Jaccard on edges), and whether the canonical regulons are recovered at all — FOXJ1's regulon must contain ciliary genes, SPDEF's must contain mucin/secretory genes. A network failing that sanity check is not used.

---

## 7. Module 6 — VIPER

- **Implementation:** `pyviper` (Califano lab Python port, for the single-cell/metacell path) and R `viper` for msVIPER on the pseudobulk contrasts; both are run and required to agree within rank tolerance. `decoupler`'s multivariate methods (mlm/ulm/wsum) are run as a third, non-VIPER scoring family for orthogonal comparison.
- **msVIPER** on each pseudobulk GES (club, goblet, ciliated vs basal), with the null from sample-label permutation (donor-level shuffling, ≥1,000 permutations) — a gene-shuffling null inflates significance and is not used for the reported statistics.
- **metaVIPER** across the per-dataset ARACNe networks and across the per-library ChEA3 networks; this is the primary ranking, since no single network is assumed correct.
- **Single-cell VIPER** on metacells to obtain a TF-activity matrix; then (a) differential activity along the basal→club pseudotime, and (b) clustering in activity space as a check that `club_core` separates from `goblet_core` on regulatory activity, not just on marker expression.
- **NaRnEA** as an alternative to VIPER's NES where regulon-size bias is a concern; reported as a secondary column.
- **Outputs per contrast:** TF, NES, p, FDR (BH), rank, regulon source, regulon size, leading-edge targets, TF expression fraction in basal / club / goblet / ciliated, and a flag for whether the TF's own transcript is detectable (VIPER's value is that it can call TFs whose mRNA is low — but that fact should be visible, not buried).

---

## 8. Module 7 — The positive-control benchmark (the gate)

`config/benchmark.yaml`, frozen before tuning:

**Gating controls** — well-powered lineages; these decide whether the pipeline runs at all.

| Contrast | Must-rank TF | Should-also-rank |
|---|---|---|
| basal → goblet | **SPDEF** | FOXA3, CREB3L1, XBP1, ELF3, EHF |
| basal → ciliated | **FOXJ1** | RFX2, RFX3, TP73, MYB, E2F4, GMNC |
| basal → deuterosomal | **MCIDAS** | FOXN4, MYB, E2F4, E2F5, GMNC |
| basal identity (reverse contrast) | **TP63** | KLF5, SOX2, PAX9, TFAP2A |

**Reported controls** — rare lineages; scored and reported, but *not* part of the pass/fail gate.

| Contrast | Must-rank TF | Should-also-rank | Power problem |
|---|---|---|---|
| basal → ionocyte | **FOXI1** | ASCL3 | Ionocytes are ~0.1–1% of airway epithelium; a single dataset rarely has enough |
| basal → tuft/brush | **POU2F3** | ASCL2, SPIB, GFI1B | Rarer still, and inconsistently captured across dissociation protocols |
| basal → PNEC | **ASCL1** | NEUROD1, INSM1 | Rare and spatially clustered, so abundance is high-variance across donors |
| basal → hillock | KLF4 | GRHL1, GRHL3, SOX15 | Human hillock identity is less settled than the mouse definition |

**Rare-lineage handling — you're right that ionocytes may not be recoverable.** Rather than dropping them, the design makes their evaluation conditional and honest:
- Each rare contrast is evaluated **only on the pooled atlas** (never per-dataset) and **only if** the lineage has ≥150 cells across ≥3 donors after annotation; otherwise the result is recorded as `underpowered`, not as a failure. The realized cell counts are reported either way, so a missing control is visibly a sampling problem rather than a silent omission.
- Where the count clears the floor, the contrast runs on metacells rather than single cells to stabilize the signature.
- Rare-lineage outcomes never enter the sweep objective. If FOXI1 recovery were allowed to influence configuration selection, a handful of ionocytes would end up steering the club analysis.
- Interpretation of a *reported* control differs from a gating one: recovering FOXI1 from 200 ionocytes is positive evidence that VIPER is working at low n; failing to is uninformative. This asymmetry is stated in the report so nobody reads a blank cell as a negative result.

The aggregate metric is a **lineage-TF recovery curve** across all gating contrasts (rank percentile of each must-rank TF, plus the fraction of should-also-rank TFs in the top 10%), with rare lineages plotted separately and marked by their cell counts.

**Acceptance criteria** (your framing — the control TF need not be #1, but must be clearly elevated):
- Primary (gating controls only): SPDEF and FOXJ1 each in the **top 5% of the tested TF universe** and in the **top 25 absolute**, with FDR < 0.05 and positive NES, in ≥2 of the 3 GES contrasts and in **both** network families. MCIDAS and TP63 in the top 10% under the same conditions.
- Secondary: ≥3 of the "should-also-rank" TFs per contrast in the top 10%.
- Specificity: SPDEF must rank *higher in the goblet contrast than in the club contrast* and FOXJ1 higher in the ciliated contrast than in either secretory contrast — this is the discriminative test that a generic "differentiation" signal would fail.
- Reported as a rank-percentile table plus recovery curves, not as a binary.

**Sweep.** A bounded grid over the hyperparameters that plausibly matter — integration method, metacell size, ARACNe p-threshold and regulon-size band, GES construction (endpoint vs trajectory), prior-network weighting, VIPER null — scored by the composite control-recovery metric. The winning configuration is written to `config/frozen.yaml`, and the club analysis runs **only** under that frozen config. The sweep results table is committed, so the selection is auditable and the degree of configuration-sensitivity is visible.

**Failure branch.** If no configuration passes, the pipeline stops and reports rather than proceeding to the club prediction. Documented diagnostic ladder: (1) is the goblet population real and correctly separated? (2) does SPDEF's regulon contain its known targets in this network? (3) is the GES dominated by ambient/secreted transcripts? (4) is the donor-pairing collapsing the contrast?

---

## 9. Module 8 — The club prediction

Run under the frozen config only:

- Ranked TF table for **club_core vs basal**, across both network families and all three GES contrasts, with a **consensus rank** (rank-aggregation across networks/contrasts, e.g. RRA or Stouffer on z-scores) and a per-TF stability score (fraction of leave-one-dataset-out runs in which the TF stays in the top 50).
- **Goblet-subtraction view:** the club-vs-basal ranking with the goblet-vs-basal ranking regressed out, plus a direct **club vs goblet** contrast. A TF that is high in both secretory contrasts is a *secretory-program* TF; a TF high in club and low/negative in goblet is a *club-specifying* candidate. Both lists are reported, labeled distinctly. This is where the annotation work pays off — without the club/goblet separation these are indistinguishable.
- **Regulator hierarchy:** for the top candidates, use ChEA-KG's signed TF–TF edges plus the ARACNe TF–TF subnetwork to order them into upstream/downstream layers, and compute VIPER's shadow/synergy analysis to remove TFs that rank highly only through overlapping regulons.
- **Chromatin corroboration — deferred to v2 (see §10a).** Not built in the first pass.
- **Deliverable table:** ranked candidates with NES, consensus rank, stability, leading-edge club targets (does the TF's leading edge contain SCGB1A1/SCGB3A2/BPIFA1?), goblet-contrast NES, expression pattern, and a short annotation of prior literature support. Explicitly separates *recovered known biology* from *novel nominations*, and flags the plausible-but-unproven ones as the experimental shortlist (CRISPRi/overexpression in ALI, readout SCGB1A1/SCGB3A2 vs MUC5AC).

Prior expectations for club biology (NKX2-1, FOXA2, SOX2, CEBP family, ELF3/EHF, XBP1/CREB3L1 secretory machinery, NOTCH effectors HES1/HEY1) are used **for interpretation only** and are deliberately not part of the tuning objective — putting them in the benchmark would make the club result circular.

---

## 9a. Deferred to v2 — other omics layers

Explicitly out of scope for the first release, but the interfaces are designed now so adding them later is not a rewrite:

- **scATAC / multiome.** A `network/atac.py` slot that produces a third regulon set (motif accessibility in club-specific peaks, SCENIC+ style) which drops into metaVIPER alongside ARACNe and the prior networks — the same interface, a different edge source. Candidates supported by both expression-inferred and accessibility-inferred regulons are the strongest nominations.
- **CUT&RUN / ChIP for the top candidates.** Turns a predicted regulon into a measured one for a handful of TFs; the converter already accepts an arbitrary signed edge list, so this enters as another `regulons.yaml` source.
- **Proteomics / secretome.** Relevant because the club phenotype is defined by what it secretes; a downstream check rather than a network input.
- **Spatial.** Would test whether the predicted club program is positionally graded along the proximo-distal axis — `GSE143868`'s 35 locations make that testable with expression alone first.

Nothing in the v1 design assumes a single omics layer: the regulon interface is `(TF, target, mode, likelihood)` regardless of where the edge came from, which is the reason multi-omics can be added without touching the VIPER, benchmark, or reporting modules.

## 10. Robustness and reporting

- Leave-one-dataset-out and leave-one-system-out (in vivo vs ALI) re-runs; report rank stability.
- Cell-subsampling bootstrap (e.g. 20× at 80%) for confidence intervals on NES.
- Null-regulon control: degree-preserving edge shuffling; the control TFs must lose significance.
- Comparison of the two network families as a headline figure — concordance of the club ranking, with the TFs unique to each called out.
- `snakemake --report` HTML + a `docs/RESULTS.md` narrative; all tables as CSV in `results/tables/`; provenance JSON (dataset checksums, network versions, config hash, package versions) attached to every output.
- CI: unit tests on the resolver rules and the regulon converter against a small synthetic fixture; a `--dry-run` DAG check; no full pipeline in CI.

---

## 11. Build order

1. Registry + loaders + gene harmonization; atlas assembled and QC'd.
2. Integration + hierarchical annotation; **secretory resolver** and its held-out validation. *(Highest-risk step — budget the most time here.)*
3. GES construction for all four contrasts.
4. Prior-network converter (ChEA3/CollecTRI → signed, weighted regulons) — deliverable on its own, reusable outside this project.
5. ARACNe runs, per-dataset and pooled.
6. VIPER/metaVIPER wrappers; benchmark harness; sweep; freeze config.
7. Club prediction, robustness, report.

**Known failure modes to watch:** ambient SCGB1A1/MUC5AC contaminating the club call; ALI cultures dominating the transition signal and biasing toward culture-stress TFs; regulon-size bias inflating broad-regulon TFs; datasets whose "club" labels are actually a merged secretory class silently defining the contrast; and the basal→secretory signature being dominated by a generic differentiation/proliferation-exit program that buries the lineage-specific TFs (the goblet-subtraction view is the mitigation).

---

## 12. Settled and open

**Settled:** repo name `Basal_to_club_TF`; dual-network design (ARACNe + ChEA3/CollecTRI prior); `GSE121600` (Ruiz García) as trajectory and resolver ground truth; HLCA core as the label-transfer reference; full lineage-TF benchmark including the rare compartments, with ionocyte/tuft/PNEC reported-but-not-gating; other omics layers deferred to v2 with the interfaces stubbed.

**Open:**
1. Confirm or edit the recommended dataset set in `docs/DATASETS.md` (the ✅ rows).
2. Held-out sample for resolver validation — default is the `GSE121600` in vivo bronchial biopsy; alternative is `madissoon_2023` for its explicit submucosal-gland labels.
3. Any in-house HBEC data to add as an anchor batch (not required).
4. Repo visibility at first push, and whether you want the GitHub Actions CI on from day one.
