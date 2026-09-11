# Basal_to_club_TF

Inference of candidate transcription factors driving the **basal → club** transition
in human airway epithelium, using VIPER/metaVIPER over an aggregated public
single-cell RNA-seq atlas.

The central problem this repository is built around is that a club-cell TF list is
easy to produce and hard to trust. Two things make it untrustworthy: the club
compartment is routinely mislabelled in published annotations (collapsed with
goblet, or split on a marker that does not separate them), and a TF ranking has no
natural null. This pipeline addresses both structurally rather than by inspection.

## The two design commitments

**1. The club compartment is re-derived, not inherited.**
Author annotations that call cells "secretory", "club" or "club/goblet" are treated
as one opinion among three. Cells are re-partitioned from expression by an explicit
resolver (`src/basal_to_club/secretory/resolver.py`) that removes submucosal-gland
and ciliated contaminants *before* the club/goblet call, never uses MUC5B as a
discriminant (it is expressed by club, goblet and gland mucous cells in human
airway), and lets a convincing goblet program veto a club call but not the reverse.
The asymmetry is deliberate: the analysis must not gain club cells by absorbing
goblet cells.

**2. The club prediction is gated on recovering known biology first.**
Before any club result is computed, the same pipeline must recover **SPDEF** for
basal→goblet and **FOXJ1** for basal→ciliated, in both network families, with a
specificity requirement that each control ranks higher in its own contrast than in
the club contrast. The thresholds live in `config/benchmark.yaml` and are frozen
before any hyperparameter sweep. If the gate fails, `rule club_prediction` does not
run — a failed gate produces a diagnostic report instead of a TF list.

```
                                          ┌── gate FAILED → benchmark report only
atlas → resolver → signatures → VIPER ────┤
                                          └── gate PASSED → club TF ranking
```

## Quickstart

```bash
make setup                 # conda env + editable install
make test                  # unit tests; needs no data
snakemake -n               # dry run of the full DAG

make atlas                 # stop after the resolved, annotated atlas
make sweep                 # grid search on the controls only -> config/frozen.yaml
make benchmark             # stop at the gate; writes the benchmark report
make all                   # full run; club prediction only if the gate passes
```

Nothing is downloaded until you run a target that needs it. Every fetch writes a
`results/provenance/*.json` with URL, sha256 and retrieval time.

## Repository layout

```
config/       datasets, markers/lineage TFs, run params, frozen benchmark, regulon sources
workflow/     Snakefile + per-stage rules
src/          the package: ingest, qc, integrate, annotate, secretory, trajectory,
              signature, network, viper, benchmark, report
tests/        logic tests that run without data
docs/         DESIGN.md (full spec), DATASETS.md (accessions and why), DECISIONS.md
```

## Data

The shortlist, with accessions, sample counts, per-study role and per-study
caveats, is in [docs/DATASETS.md](docs/DATASETS.md). The minimal recommended set is
configured in `config/datasets.yaml`.

**GSE121600** (Ruiz García et al.) is the anchor. It is an ALI differentiation time
course rather than a steady-state snapshot, so the basal→club transition is
*observed over time* rather than inferred from a static contrast. It supplies the
resolver's held-out validation sample and an independent time-course check: club
fraction must rise with ALI day. Its non-human samples are excluded in
`config/datasets.yaml` and the exclusion fails loudly if the sample identifiers
ever stop matching.

## Regulons

Two independent families, kept separate all the way through:

- **de novo** — ARACNe-AP over airway metacells built within donor.
- **prior** — ChEA3 libraries and CollecTRI.

ChEA3 libraries are unsigned gene sets; VIPER needs a mode of regulation and a
likelihood. Both are estimated from metacell TF–target correlation rather than
assumed, weak edges are **kept with likelihood 0** rather than dropped (dropping
them would quietly turn the prior network into a data-derived one), and the sign
heuristic is audited against CollecTRI's curated signs. The audit is reported, not
used to correct the network: if concordance is near chance, that is a finding about
the build and belongs in the report.

Every network passes a sanity gate before VIPER sees it — FOXJ1 must regulate cilia
genes, SPDEF must regulate mucins. A network that fails this is broken in a way no
downstream statistic makes obvious.

## Rare lineages

Ionocyte (FOXI1), tuft (POU2F3) and PNEC (ASCL1) controls are scored and reported
but **can never gate and never enter the sweep objective** — enforced in code, in
`score_gating` versus `score_reported_only`, not by convention. Below the cell and
donor floors in `config/benchmark.yaml` a contrast is recorded as `underpowered`
with its realized counts. The interpretation is asymmetric and stated in the
output: a miss is uninformative at these cell numbers, a hit is positive evidence.

## Outputs

| path | contents |
|---|---|
| `results/sweep/sweep_results.csv` | one row per configuration, scored on the controls |
| `results/atlas/atlas_resolved.h5ad` | integrated, annotated, club/goblet-resolved atlas |
| `results/secretory/resolver_holdout_metrics.json` | held-out balanced accuracy, ALI time-course trend, composition by medium |
| `results/networks/network_qc.csv` | per-network regulon sanity results |
| `results/benchmark/benchmark_report.md` | control recovery, specificity, gate verdict |
| `results/club/club_tf_ranking.csv` | ranked candidates with per-family NES, FDR, tier |
| `results/club/robustness.csv` | rank stability across resampling runs |
| `results/report.md` | final report |

## What this does not establish

VIPER infers TF activity from the expression of a TF's targets. A TF whose regulon
overlaps the club program will score highly whether or not it is causal in this
system, and prior networks are biased toward well-studied TFs. A high-confidence
candidate here is a **nomination for perturbation**, not a demonstrated driver.
Robustness columns are empty when resampling has not been run; that means stability
is unassessed, not good.

## Citation

Datasets are cited individually in `docs/DATASETS.md` and echoed into
`results/provenance/`. Methods: VIPER (Alvarez et al., Nat Genet 2016), metaVIPER
(Ding et al., Nat Commun 2018), ARACNe-AP (Lachmann et al., Bioinformatics 2016),
ChEA3 (Keenan et al., Nucleic Acids Res 2019), CollecTRI (Müller-Dott et al.,
Nucleic Acids Res 2023).

## License

MIT.

## Run outputs and QC

- [`docs/RUN_REPORT.md`](docs/RUN_REPORT.md) — what the first real run did, and the four defects it exposed
- [`docs/qc_2026-09/QC_REPORT.md`](docs/qc_2026-09/QC_REPORT.md) — QC tables and figures: study inventory,
  club/goblet resolver performance (leave-one-study-out), positive-control recovery, three-lineage driver
  comparison, and regulon evidence for the candidate basal to club drivers
- `docs/qc_2026-09/run_tables/` — the run's result tables (`results/` itself is git-ignored)

### Class-definition sensitivity (2026-09-10)

- [docs/atlas_label_arms/ARM_COMPARISON.md](docs/atlas_label_arms/ARM_COMPARISON.md) —
  re-run of VIPER with club/goblet taken from the HLCA's own `ann_finest_level` labels
  instead of the marker resolver, on the four studies with enough cells in both classes.
  Rankings correlate at rho 0.97 but 4 of the top 20 change; SPDEF recovery is *worse*
  with the atlas labels; 11 TFs are club-preferential under all three definitions; FOXA2
  is not one of them.

### Orthogonal time-course check with DREM (2026-09-11)

VIPER scores static cell states; DREM ([Ernst et al. 2007](https://doi.org/10.1038/msb4100115))
fits regulator-annotated bifurcations to an actual differentiation time course. Both arms use
the same ChEA3-derived regulons, so the comparison tests the *evidence* rather than the network.

- [docs/drem_ali/DREM_COMPARISON.md](docs/drem_ali/DREM_COMPARISON.md) — nasal ALI course
  (GSE121600, 10x, donors D246/D275, days 2-22)
- [docs/drem_bronchial/DREM_COMPARISON.md](docs/drem_bronchial/DREM_COMPARISON.md) +
  [METHODS.md](docs/drem_bronchial/METHODS.md) — bronchial ALI course (GSE233145, primary HBEC,
  Drop-seq, donors Donor_1/Donor_2, days 0-28), an independent replication in the
  tissue this repository is actually about

Both courses agree: the *ciliated* positive control is recovered at rank 1-3 of ~800 TFs in every
model (*FOXJ1*, with *RFX2/RFX3/TP73* behind it); goblet control recovery is set by sequencing
depth (nasal 10x recovers *SPDEF* at percentile 0.04, bronchial Drop-seq does not, at 0.4%
*MUC5AC*+ cells); club rankings are reproducible within donor across label arms
(rho 0.83 bronchial / 0.63 nasal) but not across donors (0.49 / 0.23); and **neither course
corroborates the static VIPER club-specific set** (median-ranked, p = 0.28 / 0.30). The club TFs
reproducible in both tissues are *FOXM1*, *MYBL2*, *TFDP1*, *PAX9*, *ZBTB7C* — three of them
cell-cycle regulators.

Scripts: `scripts/prep_drem_bronchial.py`, `scripts/run_drem_ali.py` (both sets, via
`BTC_DREM_SET`), `scripts/collect_drem_ali.py`, `scripts/compare_drem_viper.py`,
`scripts/plot_drem_cross_dataset.py`. DREM is vendored and patched in `external/STEM_DREM/`
(three null-pointer guards needed for batch mode; see `docs/drem_ali/DREM_COMPARISON.md`).
