"""Club vs goblet resolution.

The problem this solves: published airway annotations frequently collapse the
secretory compartment into one "secretory" or "club" label, or split it on a
single marker. A club prediction built on that label is a prediction about a
mixture. This module re-partitions the secretory compartment from expression,
using rules that are explicit and testable rather than a threshold chosen by eye.

Three design decisions worth stating:

1. Submucosal gland (SMG) cells are removed FIRST. Gland mucous cells are
   MUC5B-high and gland serous cells LTF/LYZ-high; both otherwise land in the
   secretory compartment and blur the club definition, especially in biopsies.
2. MUC5B is never a discriminant. In human proximal airway it is expressed by
   club cells, goblet cells and SMG mucous cells. It is carried as a descriptive
   covariate only.
3. Goblet vetoes club. A cell with a convincing goblet program is not club even
   if it also scores club-positive; the reverse is not applied. This asymmetry is
   deliberate - the analysis must not gain club cells by absorbing goblet cells,
   because that is precisely the error that would contaminate the TF prediction.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Literal

import numpy as np
import pandas as pd

CellClass = Literal["club", "goblet", "hybrid", "smg", "ciliated_contaminant",
                    "secretory_unassigned"]


#: Programs whose only role is to REMOVE cells. When one of these shows no
#: credible bimodality it means the contaminant population is absent, not that
#: it is hiding below a quantile - so the threshold goes to +inf and nothing is
#: flagged. Using a quantile fallback here would silently delete a fixed
#: fraction of the secretory compartment from every dataset.
CONTAMINANT_PROGRAMS = ("smg", "ciliated")


@dataclass
class ResolverThresholds:
    """Per-dataset thresholds. Fitted, not hardcoded - see `fit_thresholds`."""
    club: float
    goblet: float
    smg: float
    ciliated: float
    separation: dict[str, float] = field(default_factory=dict)
    fallback: dict[str, str] = field(default_factory=dict)
    method: str = "gmm"

    def to_dict(self) -> dict:
        return asdict(self)


def fit_gmm_threshold(scores: np.ndarray, min_separation: float = 2.0,
                      random_state: int = 0) -> tuple[float, float]:
    """Split a score distribution with a 2-component Gaussian mixture.

    Returns (threshold, separation), where `separation` is the component mean gap
    in pooled-SD units and is set to exactly 0.0 when the data give no evidence of
    two populations. The caller uses that as its fallback signal.

    Deciding bimodality by separation alone does not work: a 2-component mixture
    fitted to a SINGLE Gaussian still reports a mean gap of roughly 1.5-2.5 pooled
    SD, because the components tile one hump. Two conditions are therefore
    required before a threshold is trusted:

      1. the 2-component model beats the 1-component model on BIC, and
      2. the fitted mixture density is genuinely bimodal - it has an interior
         local minimum between the component means.

    Condition 2 is the definition of bimodality and is what a separation cutoff
    is a proxy for; testing it directly removes the proxy.
    """
    from sklearn.mixture import GaussianMixture

    x = np.asarray(scores, dtype=float).reshape(-1, 1)
    x = x[np.isfinite(x).ravel()]
    if x.shape[0] < 50:
        return (float(np.quantile(x, 0.75)) if x.size else 0.0), 0.0

    kw = dict(random_state=random_state, n_init=3)
    gm1 = GaussianMixture(n_components=1, **kw).fit(x)
    gm2 = GaussianMixture(n_components=2, **kw).fit(x)
    q75 = float(np.quantile(x, 0.75))
    if gm2.bic(x) >= gm1.bic(x):
        return q75, 0.0                      # one population; no threshold to find

    means = gm2.means_.ravel()
    sds = np.sqrt(gm2.covariances_.ravel())
    lo, hi = int(np.argmin(means)), int(np.argmax(means))
    pooled = float(np.sqrt((sds[lo] ** 2 + sds[hi] ** 2) / 2)) or 1e-9
    separation = float((means[hi] - means[lo]) / pooled)

    grid = np.linspace(means[lo], means[hi], 512)
    density = np.exp(gm2.score_samples(grid.reshape(-1, 1)))
    dip_at = int(np.argmin(density))
    is_bimodal = 0 < dip_at < len(grid) - 1     # interior minimum => true dip

    if not is_bimodal or separation < min_separation:
        return q75, (separation if is_bimodal else 0.0)

    return float(grid[dip_at]), separation


def fit_thresholds(scores: pd.DataFrame, min_separation: float = 2.0,
                   random_state: int = 0) -> ResolverThresholds:
    """Fit one threshold per program, within a single dataset.

    Thresholds are fitted per dataset because ambient level, chemistry and
    dissociation protocol shift the score scale; a global cut would silently
    reclassify whole datasets.
    """
    out, sep, fallback = {}, {}, {}
    for program in ("club", "goblet", "smg", "ciliated"):
        thr, s = fit_gmm_threshold(scores[program].to_numpy(), min_separation, random_state)
        sep[program] = s
        if s >= min_separation:
            out[program], fallback[program] = thr, "none"
        elif program in CONTAMINANT_PROGRAMS:
            # No contaminant population detected: flag nothing rather than
            # excluding the top quartile of an unimodal distribution.
            out[program], fallback[program] = float("inf"), "no_population_detected"
        else:
            out[program], fallback[program] = thr, "quantile_0.75"
    return ResolverThresholds(separation=sep, fallback=fallback, **out)


def classify(scores: pd.DataFrame, thr: ResolverThresholds,
             goblet_vetoes_club: bool = True, hybrid_class: bool = True) -> pd.Series:
    """Assign a class per cell from program scores and fitted thresholds.

    Order matters: contaminants are removed before the club/goblet call, so a
    ciliated or gland cell can never be counted as club.
    """
    club_pos = scores["club"] >= thr.club
    gob_pos = scores["goblet"] >= thr.goblet
    smg_pos = scores["smg"] >= thr.smg
    cil_pos = scores["ciliated"] >= thr.ciliated

    out = pd.Series("secretory_unassigned", index=scores.index, dtype=object)
    out[club_pos & ~gob_pos] = "club"
    out[gob_pos & ~club_pos] = "goblet"

    both = club_pos & gob_pos
    if hybrid_class:
        out[both] = "hybrid"
    elif goblet_vetoes_club:
        out[both] = "goblet"          # asymmetric on purpose: never absorb goblet into club
    else:
        out[both] = "club"

    out[smg_pos] = "smg"              # contaminant classes override
    out[cil_pos] = "ciliated_contaminant"
    return out.astype("category")


def score_programs(adata, markers: dict[str, list[str]], alias_map=None,
                   score_method: str = "score_genes") -> pd.DataFrame:
    """Score the four programs the resolver needs, on log-normalized data.

    SMG is scored as the max of the serous and mucous panels: either one is
    sufficient grounds for exclusion, and averaging them would let a strongly
    serous cell slip through.
    """
    import scanpy as sc

    from basal_to_club.utils.genes import require_panel

    panels = {
        "club": markers["club"],
        "goblet": markers["goblet"],
        "ciliated": markers["ciliated"],
        "smg_serous": markers["smg_serous"],
        "smg_mucous": markers["smg_mucous"],
    }
    scores = {}
    for name, panel in panels.items():
        genes = require_panel(name, [g for g in panel if g != "MUC5B"],
                              adata.var_names, alias_map=alias_map)
        sc.tl.score_genes(adata, genes, score_name=f"_score_{name}", use_raw=False)
        scores[name] = adata.obs.pop(f"_score_{name}")

    df = pd.DataFrame(scores, index=adata.obs_names)
    df["smg"] = df[["smg_serous", "smg_mucous"]].max(axis=1)
    return df[["club", "goblet", "ciliated", "smg"]]


def resolve_secretory(adata, markers: dict, params: dict, dataset_key: str = "dataset_id"):
    """Full resolver: score -> fit per dataset -> classify. Returns (assignments, thresholds)."""
    sec = params["secretory"]
    scores = score_programs(adata, markers, score_method=sec.get("score_method", "score_genes"))
    assignments, fitted = [], {}
    for ds, idx in adata.obs.groupby(dataset_key, observed=True).groups.items():
        sub = scores.loc[idx]
        thr = fit_thresholds(sub, sec.get("min_component_separation", 2.0))
        fitted[str(ds)] = thr.to_dict()
        assignments.append(classify(sub, thr,
                                    goblet_vetoes_club=sec.get("goblet_vetoes_club", True),
                                    hybrid_class=sec.get("hybrid_class", True)))
    out = pd.concat(assignments).reindex(adata.obs_names)
    return pd.DataFrame({"cell_class": out}).join(scores), fitted


def balanced_accuracy(truth: pd.Series, pred: pd.Series,
                      classes: tuple[str, ...] = ("club", "goblet")) -> float:
    """Balanced accuracy on the holdout, restricted to the classes under test."""
    from sklearn.metrics import balanced_accuracy_score

    mask = truth.isin(classes) & pred.isin(classes)
    if mask.sum() == 0:
        return float("nan")
    return float(balanced_accuracy_score(truth[mask], pred[mask]))


def monotonic_club_fraction(obs: pd.DataFrame, day_key: str = "ali_day",
                            class_key: str = "cell_class") -> dict:
    """Independent check from the GSE121600 ALI time course.

    Club abundance should rise with ALI day. A resolver that reproduces holdout
    labels but shows no time trend is splitting a continuum at an arbitrary point.
    """
    from scipy.stats import spearmanr

    sub = obs.dropna(subset=[day_key])
    if sub.empty:
        return {"status": "no_timecourse_samples"}
    frac = (sub.groupby(day_key, observed=True)[class_key]
               .apply(lambda s: float((s == "club").mean())))
    if frac.size < 3:
        return {"status": "too_few_timepoints", "n_timepoints": int(frac.size)}
    rho, p = spearmanr(frac.index.astype(float), frac.to_numpy())
    return {"status": "ok", "spearman_rho": float(rho), "p_value": float(p),
            "club_fraction_by_day": {float(k): float(v) for k, v in frac.items()},
            "monotonic_increase": bool(rho > 0 and p < 0.05)}
