# Decision log

Choices that a reviewer would otherwise have to reverse-engineer from the code,
with the reasoning and the alternative that was rejected.

## D1 — The club compartment is re-derived rather than inherited

*Decision.* Author cell-type labels enter as one of three opinions (author,
reference-transferred, marker-based), arbitrated by modal vote; the secretory
compartment is then re-partitioned by a dedicated resolver.

*Why.* Published airway annotations disagree about what "club" means. Some collapse
club and goblet into "secretory", some split them on a single marker. A TF
prediction built on an inherited label is a prediction about whatever mixture that
label denotes.

*Rejected.* Trusting the finest available label. It is not comparable across the
studies being aggregated, which is the reason for aggregating in the first place.

## D2 — MUC5B is never a club/goblet discriminant

*Decision.* MUC5B is carried as a descriptive covariate only. A unit test enforces
that it cannot separate the two panels.

*Why.* In human proximal airway, MUC5B is expressed by club cells, goblet cells and
submucosal-gland mucous cells. Using it as a discriminant would create a split that
looks clean and is not.

## D3 — Gland and ciliated cells are removed before the club/goblet call

*Decision.* SMG (serous or mucous) and ciliated contaminants are classified first
and override the secretory call. SMG is scored as the **max** of the serous and
mucous panels, not their average.

*Why.* Gland cells otherwise land in the secretory compartment and blur the club
definition, particularly in biopsy datasets. Averaging the two gland panels would
let a strongly serous cell slip through.

## D4 — Goblet vetoes club; club does not veto goblet

*Decision.* When both programs are positive the cell is called `hybrid`; if hybrids
are disabled, it is called goblet.

*Why.* The failure mode that would corrupt the result is club gaining goblet cells,
because the club TF prediction would then be partly a goblet prediction. The
reverse error is much less costly here, so the rule is asymmetric on purpose.

## D5 — An absent contaminant population must not be invented

*Decision.* Thresholds are fitted per dataset by 2-component Gaussian mixture, and
bimodality is tested by (a) BIC preference for two components over one and (b) an
interior local minimum in the fitted density. For the contaminant channels (SMG,
ciliated), failing that test sets the threshold to `+inf` — nothing is flagged.

*Why.* Two separate errors were possible here, and the first version had both. A
2-component mixture fitted to a single Gaussian still reports a mean gap of roughly
1.5–2.5 pooled SD, so a separation cutoff does not distinguish "two populations"
from "one hump tiled by two components". And a quantile fallback on a contaminant
channel deletes a fixed fraction of the secretory compartment from every dataset
whether or not gland cells are present. Both are covered by tests
(`test_unimodal_scores_do_not_get_a_confident_split`,
`test_absent_contaminant_population_does_not_delete_cells`).

## D6 — Batch-corrected expression never reaches VIPER

*Decision.* The integrated latent space is used for neighbours, clustering and
label transfer only. Every signature is recomputed from raw counts, pseudobulked
within donor.

*Why.* Batch correction moves genes in ways a regulon-based method reads as
regulatory signal. The correction is needed to define the cells, not to compare
them.

## D7 — Donor is the unit of replication

*Decision.* Signatures come from donor-level pseudobulk with limma-voom and donor
blocking, not from cell-level differential tests.

*Why.* Cell-level tests treat thousands of cells from one donor as independent; the
resulting p-values are largely a function of cell count, and VIPER inherits that
inflation directly.

## D8 — Two network families, kept separate

*Decision.* De novo (ARACNe) and prior-knowledge (ChEA3, CollecTRI) regulons are
run independently, and a gating control must pass in **both**.

*Why.* Agreement between an expression-inferred network and a curated one is the
evidence that a result is not an artifact of one network construction. Merging them
into a single edge list would destroy exactly that comparison.

## D9 — Weak prior edges are kept with likelihood 0

*Decision.* Prior edges whose TF–target correlation falls below `min_abs_rho` stay
in the regulon with likelihood 0 rather than being removed.

*Why.* Removing them makes regulon membership itself data-dependent, which turns
the prior network into a second de novo network and collapses D8.

## D10 — Inferred signs are audited, not corrected

*Decision.* Signs assigned to ChEA3 edges are compared with CollecTRI's curated
signs and the concordance is reported.

*Why.* If concordance is near chance the heuristic is not working on this dataset.
That is a finding about the build; silently overwriting signs would hide it.

## D11 — The benchmark gates execution, and its thresholds are frozen first

*Decision.* `config/benchmark.yaml` is written before any sweep. The sweep is scored
on `gating_recovery_score`, which is computed from the gating contrasts only. The
winning configuration is frozen to `config/frozen.yaml`, and only then is the club
contrast evaluated.

*Why.* Selecting a configuration on a metric that can see the club result is
circular. Freezing first is what makes "the controls came out where they should" a
claim rather than a consequence of tuning.

The sweep itself is scoped to match: the grid contains only parameters downstream
of integration, `run_sweep` takes no club-related input, and `sweep.py` never
computes the club or rare-lineage contrasts. It runs on one representative network
per family — the sweep compares configurations, not networks, and the gate itself
still requires both families at full network depth.

*Rejected.* Relaxing thresholds when the gate fails. `check_gate.py` says so
explicitly in its failure message.

## D12 — Specificity is part of the gate

*Decision.* Each control TF must rank higher in its own contrast than in the club
contrast.

*Why.* Without it, a pipeline that merely detects "no longer basal" passes the
benchmark and then returns generic differentiation factors as club drivers.

## D13 — Rare lineages report, never gate

*Decision.* Ionocyte, tuft and PNEC controls are scored by the same code but
returned in a separate frame that the sweep objective cannot see. Below the cell and
donor floors they are recorded as `underpowered` with realized counts.

*Why.* A handful of ionocytes cannot support a stable NES, and letting such a
contrast influence configuration selection would inject noise into every downstream
choice. Recording the miss with its counts keeps a blank cell from being read as a
negative result.

## D14 — Other omics layers are deferred, but the interface is not

*Decision.* No scATAC/multiome, CUT&RUN or proteomics in v1. The regulon interface
is `(tf, target, mo, likelihood)` regardless of edge provenance.

*Why.* A chromatin-derived network can be added later as another regulon file
without touching the VIPER, benchmark or reporting modules.

## D15 — GSE121600 as anchor, with medium as a covariate

*Decision.* The Ruiz García ALI time course supplies the held-out resolver
validation sample and an independent monotonicity check. Its two culture media are
treated as a covariate with a sensitivity run, not as noise to average over.

*Why.* It observes differentiation in real time rather than inferring it from a
steady-state snapshot, which is what the transition of interest requires. Culture
medium measurably shifts secretory composition, so pooling the media without
checking would confound the very compartment being resolved.
