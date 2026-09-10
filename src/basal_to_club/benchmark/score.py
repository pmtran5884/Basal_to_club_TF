"""Scoring of the frozen positive-control benchmark.

The thresholds live in config/benchmark.yaml and are frozen before any sweep, so
this module only reads them. Two properties are enforced structurally rather than
by convention:

1. Rare-lineage contrasts are scored by the same code but returned in a separate
   frame, and `gating_recovery_score` cannot see them. A handful of ionocytes can
   therefore never steer configuration selection.
2. The specificity test is part of the gate. A TF that ranks high in every
   contrast is detecting "differentiation", not "goblet". SPDEF must rank higher
   in basal->goblet than in basal->club for the goblet control to count as passed.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd


@dataclass
class ControlResult:
    contrast: str
    tf: str
    kind: str                 # "primary" | "secondary"
    rank: int | None
    rank_percentile: float | None
    nes: float | None
    fdr: float | None
    n_tf_universe: int
    passed: bool
    reason: str

    def to_dict(self) -> dict:
        return asdict(self)


def rank_tf(activity: pd.DataFrame, contrast: str, tf: str) -> dict:
    """Rank a TF within one contrast. Rank 1 = highest positive activity.

    Ranking is on NES descending, so a TF with strongly NEGATIVE activity ranks
    last rather than first: the controls are activators of their lineage, and a
    repressed TF is not evidence for the pipeline working.
    """
    sub = activity[activity["contrast"] == contrast]
    if sub.empty:
        return {"rank": None, "rank_percentile": None, "nes": None, "fdr": None, "n": 0}
    ordered = sub.sort_values("nes", ascending=False).reset_index(drop=True)
    hit = ordered.index[ordered["tf"] == tf]
    n = len(ordered)
    if len(hit) == 0:
        return {"rank": None, "rank_percentile": None, "nes": None, "fdr": None, "n": n}
    i = int(hit[0])
    row = ordered.loc[i]
    return {"rank": i + 1, "rank_percentile": (i + 1) / n,
            "nes": float(row["nes"]), "fdr": float(row.get("fdr", np.nan)), "n": n}


def score_contrast(activity: pd.DataFrame, contrast: str, spec: dict,
                   requirements: dict) -> list[ControlResult]:
    """Score the primary and secondary TFs of one contrast."""
    out: list[ControlResult] = []
    max_pct = spec.get("max_rank_percentile", 0.05)
    max_abs = spec.get("max_absolute_rank", 25)
    max_fdr = requirements.get("max_fdr", 0.05)
    need_pos = requirements.get("require_positive_nes", True)

    primary = spec["primary"]
    r = rank_tf(activity, contrast, primary)
    reasons = []
    if r["rank"] is None:
        passed, reasons = False, ["TF absent from the tested universe"]
    else:
        if r["rank_percentile"] > max_pct:
            reasons.append(f"rank percentile {r['rank_percentile']:.3f} > {max_pct}")
        if r["rank"] > max_abs:
            reasons.append(f"absolute rank {r['rank']} > {max_abs}")
        if need_pos and not (r["nes"] > 0):
            reasons.append(f"NES {r['nes']:.2f} not positive")
        if np.isfinite(r["fdr"]) and r["fdr"] > max_fdr:
            reasons.append(f"FDR {r['fdr']:.3g} > {max_fdr}")
        passed = not reasons
    out.append(ControlResult(contrast, primary, "primary", r["rank"], r["rank_percentile"],
                             r["nes"], r["fdr"], r["n"], passed,
                             "; ".join(reasons) or "ok"))

    for tf in spec.get("secondary", []) or []:
        s = rank_tf(activity, contrast, tf)
        in_top = (s["rank_percentile"] is not None and s["rank_percentile"] <= 0.10)
        out.append(ControlResult(contrast, tf, "secondary", s["rank"], s["rank_percentile"],
                                 s["nes"], s["fdr"], s["n"], bool(in_top),
                                 "in top 10%" if in_top else "not in top 10%"))
    return out


def check_specificity(activity: pd.DataFrame, rules: dict) -> pd.DataFrame:
    """A control TF must rank higher in its OWN contrast than in the club contrast.

    Without this a pipeline that simply detects "not basal any more" would pass
    the benchmark and then hand back a club TF list made of generic
    differentiation factors.
    """
    rows = []
    for tf, rule in (rules or {}).items():
        own = rank_tf(activity, rule["must_rank_higher_in"], tf)
        club = rank_tf(activity, rule["than_in"], tf)
        if own["rank"] is None or club["rank"] is None:
            rows.append({"tf": tf, "own_contrast": rule["must_rank_higher_in"],
                         "own_rank": own["rank"], "club_rank": club["rank"],
                         "passed": False, "reason": "TF missing from one contrast"})
            continue
        ok = own["rank"] < club["rank"]
        rows.append({"tf": tf, "own_contrast": rule["must_rank_higher_in"],
                     "own_rank": own["rank"], "club_rank": club["rank"], "passed": bool(ok),
                     "reason": "specific" if ok else
                               "ranks no better in its own contrast than in club - "
                               "signal is generic differentiation, not lineage identity"})
    return pd.DataFrame(rows)


def score_gating(activity_by_family: dict[str, pd.DataFrame], bench: dict) -> tuple[pd.DataFrame, dict]:
    """Score every gating contrast in every network family, then apply the gate.

    `activity_by_family` maps family name ("de_novo", "prior") -> activity table.
    A primary control must pass in BOTH families: agreement between an
    expression-inferred network and a prior-knowledge network is the evidence that
    the result is not an artifact of one network construction.
    """
    gating = bench["gating"]
    req = gating["requirements"]
    frames = []
    for family, act in activity_by_family.items():
        for contrast, spec in gating["contrasts"].items():
            for res in score_contrast(act, contrast, spec, req):
                row = res.to_dict()
                row["family"] = family
                frames.append(row)
    scores = pd.DataFrame(frames)

    primary = scores[scores["kind"] == "primary"]
    per_contrast = {}
    for contrast in gating["contrasts"]:
        sub = primary[primary["contrast"] == contrast]
        families_passed = sorted(sub.loc[sub["passed"], "family"].unique())
        need_both = req.get("require_both_network_families", True)
        ok = (len(families_passed) == len(activity_by_family)) if need_both else bool(families_passed)

        sec = scores[(scores["contrast"] == contrast) & (scores["kind"] == "secondary")]
        sec_frac = float(sec["passed"].mean()) if len(sec) else 1.0
        if sec_frac < req.get("min_secondary_frac_in_top10pct", 0.5):
            ok = False
        per_contrast[contrast] = {
            "passed": bool(ok), "families_passed": families_passed,
            "secondary_frac_in_top10pct": sec_frac,
            "failure_reasons": sorted(set(sub.loc[~sub["passed"], "reason"])) if not ok else [],
        }

    spec_df = pd.concat(
        [check_specificity(act, req.get("specificity", {})).assign(family=fam)
         for fam, act in activity_by_family.items()], ignore_index=True) \
        if req.get("specificity") else pd.DataFrame(columns=["passed"])
    spec_ok = bool(spec_df["passed"].all()) if len(spec_df) else True

    gate = {
        "passed": bool(all(v["passed"] for v in per_contrast.values()) and spec_ok),
        "per_contrast": per_contrast,
        "specificity_passed": spec_ok,
        "specificity_detail": spec_df.to_dict("records") if len(spec_df) else [],
        "gating_recovery_score": gating_recovery_score(scores),
    }
    return scores, gate


def gating_recovery_score(scores: pd.DataFrame) -> float:
    """Single scalar used as the sweep objective. Gating contrasts only.

    Mean over primaries of (1 - rank_percentile), plus a bonus for secondaries in
    the top 10%. Rare lineages are structurally absent: they are never in the
    frame this function receives.
    """
    prim = scores[scores["kind"] == "primary"].dropna(subset=["rank_percentile"])
    if prim.empty:
        return 0.0
    base = float((1.0 - prim["rank_percentile"]).mean())
    sec = scores[scores["kind"] == "secondary"]
    bonus = 0.1 * float(sec["passed"].mean()) if len(sec) else 0.0
    return round(base + bonus, 6)


def score_reported_only(activity: pd.DataFrame, bench: dict,
                        cell_counts: dict[str, dict]) -> pd.DataFrame:
    """Score rare-lineage contrasts. Never gates, never enters the sweep objective.

    Each contrast is evaluated only if the lineage clears the cell and donor
    floor. Otherwise the row is recorded as `underpowered` with the realized
    counts, so a missing control is visibly a sampling problem and not a silent
    omission. A miss here is uninformative; a hit is positive evidence that VIPER
    works at low n. That asymmetry is written into the output so nobody reads a
    blank cell as a negative result.
    """
    cfg = bench["reported_only"]
    rows = []
    for contrast, spec in cfg["contrasts"].items():
        lineage = contrast.replace("basal_to_", "")
        counts = cell_counts.get(lineage, {"n_cells": 0, "n_donors": 0})
        powered = (counts["n_cells"] >= cfg["min_cells"]
                   and counts["n_donors"] >= cfg["min_donors"])
        if not powered:
            rows.append({"contrast": contrast, "tf": spec["primary"], "kind": "primary",
                         "status": "underpowered", "rank": None, "rank_percentile": None,
                         "nes": None, "n_cells": counts["n_cells"],
                         "n_donors": counts["n_donors"], "gates": False,
                         "interpretation": "not evaluated - insufficient cells/donors"})
            continue
        r = rank_tf(activity, contrast, spec["primary"])
        hit = r["rank_percentile"] is not None and r["rank_percentile"] <= 0.10
        rows.append({"contrast": contrast, "tf": spec["primary"], "kind": "primary",
                     "status": "evaluated", "rank": r["rank"],
                     "rank_percentile": r["rank_percentile"], "nes": r["nes"],
                     "n_cells": counts["n_cells"], "n_donors": counts["n_donors"],
                     "gates": False,
                     "interpretation": ("recovered - positive evidence VIPER works at low n"
                                        if hit else
                                        "not recovered - uninformative at this cell count")})
    return pd.DataFrame(rows)
