# Grant figures: the available data cannot identify the basal-to-club driver

Three figures and their legends. Every number is computed by a script in `scripts/`
from tables in `docs/grant_figures/tables/`; the provenance line under each legend
names the script and the table. Rendered at 300 dpi PNG plus vector PDF.

The argument runs: (1) no public dataset has both club depth and a differentiation
time axis; (2) the club compartment is not even consistently *defined* across labeling
methods; (3) at the cell numbers and donor counts that exist, a club TF ranking is not
reproducible, while the same pipeline recovers the ciliated positive control easily;
(4) the two independent inference methods agree on the ciliated answer and disagree
completely on club. Therefore the limiting resource is a purpose-built bronchial ALI
time course, not a better algorithm.

---

## Figure 1 — No available dataset supports the inference

**a. Club depth versus time coverage.** Each point is one dataset in the aggregated
corpus: 10 Human Lung Cell Atlas (HLCA) tissue studies (open circles) and the two
air-liquid-interface (ALI) differentiation time courses (filled). *x*, club cells in the
dataset (log scale); *y*, number of distinct differentiation time points sampled; marker
area scales with total epithelial cells. The shaded quadrant marks the minimum joint
requirement for trajectory-based TF inference used throughout this proposal: >=3,000
club cells **and** >=3 time points. No dataset falls in it. Tissue atlases have no
differentiation axis by construction (*y* = 0); the two ALI courses have a time axis
but only 2 donors each, and the deepest (GSE233145, 9,597 club cells across 7 days) is
Drop-seq at 418 median genes per cell.

**b. The same cells, three club definitions.** For each HLCA study, cells called club by
(i) the depositors' own annotation, (ii) HLCA label transfer as implemented in this
pipeline, and (iii) the marker-based club/goblet resolver. Bars are counts on a log axis;
grey annotations give the max/min ratio per study where it exceeds 5-fold. The three
definitions differ by up to 25-fold on the same cells, i.e. before any TF is inferred the
denominator of the analysis is not agreed.

**c. Club is the scarce arm in every study.** Distribution of cells per donor-sample
entering a differential contrast, club (red) versus the ciliated positive control (grey).
Solid line, the pseudobulk floor (25 cells) below which a donor-sample is dropped; dashed
line, the club median (88 cells per donor-sample). Club donor-samples are both fewer
(n = 46 versus 76) and smaller.

**d. Time-course coverage.** Donor x day grid for the two ALI datasets; point area and the
adjacent number give club cells recovered at that donor-day. GSE233145 samples the same
7 days in both donors; GSE121600 samples 8 days but only 2 of them in both donors, so
most nasal time points are single-donor observations and day cannot be separated from
donor there.

*Provenance: `scripts/plot_grant_fig1.py`; tables `dataset_landscape.csv`,
`club_counts_by_method.csv`, `ali_timecourse_coverage.csv`, `run_tables/pseudobulk_meta.csv`.*

### Methods and calculations behind Figure 1

- **Corpus.** 10 HLCA core studies plus GSE121600 (nasal ALI) and GSE233145 (bronchial
  ALI), harmonized to a common gene space and a 7-class epithelial vocabulary
  (basal, club, goblet, ciliated, suprabasal, ionocyte/tuft, other).
- **Club definition (i)**, published annotation: the depositor's `cell_type` mapped to the
  common vocabulary by a fixed synonym table.
- **Club definition (ii)**, label transfer: a logistic classifier trained on HLCA core
  epithelium (held-out accuracy 0.881) applied to each dataset's scaled expression.
- **Club definition (iii)**, resolver: a two-family marker gate requiring club-positive
  (*SCGB1A1*, *SCGB3A1*, *BPIFB1*, *TFF3*) above a per-dataset threshold **and**
  goblet-negative (*MUC5AC*, *MUC5B*, *SPDEF*); cells positive for both families are held
  out as hybrid and cells positive for neither as secretory-unassigned, so resolver counts
  are conservative by design.
- **Cells per donor-sample (panel c)** is the count actually entering `make_pseudobulk`
  after the 25-cell floor, not the raw label count.
- **Requirement quadrant (panel a)** is the empirical threshold from Figure 3a
  (>=1,500 club cells for 19/20 top-20 stability, rounded up to 3,000 to cover
  two independent halves) crossed with >=3 time points, the minimum for a monotone
  trajectory fit.

---

## Figure 2 — The pipeline works on controls and fails on club

**a. Positive-control recovery by method.** Rank percentile (0 = top of the ranking) of
canonical lineage regulators in the contrast for their own lineage, for static VIPER
(blue) and for the DREM time-course models in each dataset (nasal red, bronchial orange).
Shading marks the top 5%. Ciliated regulators (*FOXJ1*, *RFX2*, *RFX3*, *TP73*, *MYB*)
land in the top 5% in every method and both tissues — *FOXJ1* at percentile 0.001
(VIPER), 0.002 (nasal DREM), 0.003 (bronchial DREM). Goblet regulators are recovered by
VIPER (*SPDEF* 0.033, *CREB3L1* 0.017, *XBP1* 0.003) and by the nasal time course
(*SPDEF* 0.041) but not by the bronchial one (*SPDEF* 0.426), which is a depth failure:
*MUC5AC* is detected in 0.24% of bronchial Drop-seq cells versus 7.7% of nasal 10x cells.
Club priors (*ELF3*, *KLF5*, *FOXA2*, *CEBPB*, *SOX2*) are mid-pack in the time courses
(0.03-0.57) with no consistent ordering. The controls therefore certify the machinery
and localize the failure to the club arm.

**b. Cross-method agreement on the club ranking.** Spearman rho (colour) and shared
members of the top 20 (text, x/20) between the four club TF rankings the project
produced: static VIPER under each of the two labeling arms, and the DREM time course in
each tissue. The two static arms agree (rho = 0.97, 16/20 shared). The two time courses
agree weakly (rho = 0.34, 2/20). Static versus time course is indistinguishable from
zero (rho = -0.01 and 0.04; **0/20 shared in all four cross-family comparisons**).

**c. Top-10 candidate club drivers by method.** The four ranked lists side by side. No TF
appears in both a static and a time-course top-10. The static lists are dominated by
secretory/ER-stress regulators (*CREB3L4*, *XBP1*, *CREB3L2*, *FOXA1/2*), several of
which are goblet rather than club programs; the bronchial time-course list is topped by
ciliated regulators (*FOXJ1*, *CCDC17*, *RFX3*), i.e. at this depth the bronchial club
calls still carry transitional-ciliated signal.

**d. Formal convergence test.** Median rank percentile of the VIPER club TF set (blue)
versus all other TFs (grey) within the DREM club models of each dataset; 0.5 is chance.
Neither dataset shows enrichment (nasal p = 0.30, n = 8 testable TFs; bronchial p = 0.28,
n = 6; one-sided Mann-Whitney). Applying the pre-registered club-specificity rule — a TF
must appear in the club models and in neither the goblet nor the ciliated models —
returns **zero** TFs in both datasets.

*Provenance: `scripts/plot_grant_fig2.py`; tables `control_recovery_by_method.csv`,
`method_rank_agreement_spearman.csv`, `method_rank_agreement_top20.csv`,
`club_top10_by_method.csv`; `results/drem_*/drem_viper_convergence.json`.*

### Methods and calculations behind Figure 2

- **Static arm (VIPER).** Pseudobulk per donor x class; limma-voom club-versus-basal
  contrast blocked on donor; the moderated *t* statistics rank-transformed to normal
  quantiles; aREA (`pyviper`) against a regulon set built from Ma'ayan-lab TF-target
  resources intersected with the expressed gene space, with regulons restricted to
  25-500 targets. Run independently under each of the two labeling arms.
- **Time-course arm (DREM).** Pseudobulk expression series per day (log2 fold change
  versus day 0), one series per lineage x labeling arm x regulon family, fitted with the
  Dynamic Regulatory Events Miner; a TF's score is its best split-assignment strength
  across the fitted path model. 32 models per dataset (2 labeling arms x 4 lineages x
  2 regulon families x 2 depth settings).
- **Rank percentile** = rank of the TF divided by the number of TFs scored in that model,
  so values are comparable across models with different TF counts. For the time courses a
  TF's percentile is the best (smallest) value across the model arms of that dataset.
- **Agreement metrics (panel b).** Spearman rho over the intersection of TFs scored by
  both methods; top-20 overlap is the size of the intersection of the two top-20 sets.
- **Convergence test (panel d).** One-sided Mann-Whitney U of the club TF set's rank
  percentiles against the background of all other scored TFs, i.e. "does the static club
  set sit higher in the time-course club models than chance".
- **Club-specificity rule.** Present in the top decile of the club models and absent from
  the top decile of the goblet and ciliated models of the same dataset.

---

## Figure 3 — What limits the inference, and what the new experiment must deliver

**a. Rank accuracy versus available cells.** Split-half resampling of the nasal ALI course:
the basal arm is held at full size and the case arm is subsampled to a target cell count,
allocated across donor-day groups in proportion to the real composition, then pushed
through the unchanged production path (pseudobulk with 25-cell floor -> limma-voom
donor-blocked contrast -> rank normalization -> aREA). *y*, fraction of the full-coverage
top-20 TFs recovered. Club (red) reaches 19/20 only at ~1,500 case cells; the ciliated
control (blue) is already at 17/20 by 200 cells. Below 100 club cells no contrast can be
formed at all, because no donor-day group clears the pseudobulk floor. The median public
dataset offers 88 club cells per donor (dotted line).

**b. Where the variance actually is.** Spearman rho between two club TF rankings that
differ in exactly one respect: technical resampling only (split-half at fixed donor and
day composition, rho = 1.00); the cell-labeling method (rho = 0.63 nasal, 0.83
bronchial); the donor (rho = 0.23 nasal, 0.49 bronchial); the tissue (rho = 0.37); the
inference method (rho = 0.02). Sampling noise is not the bottleneck — donor identity,
the club definition, and the choice of inference method each move the ranking far more
than resampling does. This is why more cells from the same two donors would not settle
the question.

**c. Donors required.** Spearman-Brown projection of the reliability of a *k*-donor
consensus club ranking from the observed single-pair inter-donor agreement,
r_kk = k*r / (1 + (k-1)*r), evaluated at the bronchial (r = 0.493) and nasal (r = 0.228)
observed values. Reaching reliability 0.8 needs 5 donors under the optimistic estimate
and 14 under the pessimistic one. Every existing ALI time course has 2 (shaded).

**d. Specification for the proposed experiment.** Design parameters, the best value
available in public data, and the value required by panels a-c plus Figure 2a.

*Provenance: `scripts/plot_grant_fig3.py` and `scripts/power_club_subsample.py`;
tables `power_curve_summary.csv`, `agreement_ladder.csv`,
`donor_reliability_spearman_brown.csv`; `results/drem_bronchial/platform_depth_control.csv`.*

### Methods and calculations behind Figure 3

- **Subsampling design.** For each target cell count n, cells of the case lineage are
  drawn without replacement twice (disjoint halves) with per-group allocation proportional
  to the observed donor-day composition, so both halves retain the experiment's structure
  and only depth varies. 3 seeds per point. The basal reference arm is never subsampled,
  so the curve isolates case-arm depth.
- **Two distinct readouts.** *Split-half rho* (agreement of the two halves with each other)
  stays high even at small n, because both halves inherit the same donor-day bias.
  *Concordance with the full-coverage ranking* — plotted here — is the one that degrades,
  and the mechanism is group dropout: at small n only one or two donor-day groups clear
  the 25-cell floor, so the contrast describes those groups rather than the lineage.
  Reporting only split-half agreement would overstate the power of a small dataset.
- **Agreement ladder (panel b).** Each rung is a rho between two rankings produced by the
  identical pipeline, differing in one factor. The labeling, donor and tissue rungs reuse
  the reproducibility and cross-dataset records already reported for the time-course arm;
  the method rung is Figure 2b's static-versus-time-course cell.
- **Spearman-Brown (panel c).** Standard reliability projection treating a donor as a
  rater of the TF ranking; r is the observed between-donor rho at fixed labeling arm.
  It answers "how many donors before the consensus ranking is stable", and it is an
  upper bound on optimism because it assumes donors are exchangeable.
- **Depth requirement (panel d).** Set by the goblet-exclusion failure in Figure 2a:
  a 418-genes-per-cell library detects *MUC5AC* in 0.24% of cells, which is too sparse to
  separate club from goblet; the nasal 10x library at 3,574 genes per cell detects it in
  7.7%. >=3,000 median genes per cell and >=5% *MUC5AC*-positive cells are the operating
  points at which the resolver's goblet-negative gate was demonstrated to work.
- **Not assessed.** Batch structure of the proposed experiment, cost, and any
  alternative-design comparison (e.g. more time points versus more donors at fixed
  total cells) are outside these figures.
