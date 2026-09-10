"""The club TF prediction. Runs only downstream of a passed gate.

Every candidate is reported with the evidence a reader needs to judge it: the
NES in each network family, how many of the three signature constructions
support it, and whether it is one of the prior-plausible club TFs (annotation
only - that list never influenced the ranking).
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from basal_to_club.utils.io import load_yaml, setup_logging

CLUB = "basal_to_club"


def rank_club_tfs(combined: pd.DataFrame, by_family: pd.DataFrame,
                  plausible: list[str]) -> pd.DataFrame:
    club = combined[combined["contrast"] == CLUB]
    if club.empty:
        raise ValueError("no basal_to_club activity - the contrast was never computed")

    wide = club.pivot_table(index="tf", columns="method", values="nes")
    fam = (by_family[by_family["contrast"] == CLUB]
           .pivot_table(index="tf", columns="family", values="nes"))
    fdr = club.groupby("tf")["fdr"].min()

    out = wide.join(fam, rsuffix="_family").join(fdr.rename("min_fdr"))
    out["mean_nes"] = wide.mean(axis=1)
    out["n_methods_supporting"] = (wide > 0).sum(axis=1)
    out["method_consistency"] = 1 - wide.std(axis=1) / wide.abs().mean(axis=1).replace(0, np.nan)
    out["both_families_positive"] = (fam > 0).all(axis=1) if fam.shape[1] > 1 else False
    out["prior_plausible_club_tf"] = out.index.str.upper().isin(
        [p.upper().replace("_", "-") for p in plausible])

    out = out.sort_values("mean_nes", ascending=False)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    out["tier"] = np.where(
        (out["n_methods_supporting"] >= 2) & out["both_families_positive"] & (out["min_fdr"] < 0.05),
        "high_confidence",
        np.where(out["min_fdr"] < 0.1, "candidate", "weak"))
    return out


def main(sm):
    log = setup_logging(sm.log[0])
    markers = load_yaml("config/markers.yaml")
    combined = pd.read_csv(sm.input.activity)
    by_family = pd.read_csv(sm.input.by_family)

    ranking = rank_club_tfs(combined, by_family, markers.get("club_prior_plausible", []))
    ranking.to_csv(sm.output.ranking)
    n_high = int((ranking["tier"] == "high_confidence").sum())
    log.info("club prediction: %d TFs ranked, %d high-confidence", len(ranking), n_high)
    log.info("top 15: %s", ranking.head(15).index.tolist())

    top = ranking.head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, 8))
    colors = ["#B22222" if p else "#4C72B0" for p in top["prior_plausible_club_tf"]]
    ax.barh(top.index, top["mean_nes"], color=colors)
    ax.set_xlabel("mean NES, basal -> club")
    ax.set_title("Top candidate club drivers\n(red = prior-plausible; annotation only)")
    fig.savefig(sm.output.fig, dpi=200, bbox_inches="tight")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
