"""Robustness of the club ranking.

Leave-one-dataset-out, leave-one-system-out (in vivo vs ALI), bootstrap over
donors, and a shuffled-regulon null. A candidate that survives all four is
reportable; one that appears only with a particular dataset present is a property
of that dataset and is labelled as such rather than dropped silently.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from basal_to_club.utils.io import setup_logging


def summarize_stability(rank_matrix: pd.DataFrame, top_n: int = 25) -> pd.DataFrame:
    """rank_matrix: TFs x resampling runs, values are ranks."""
    in_top = (rank_matrix <= top_n)
    return pd.DataFrame({
        "median_rank": rank_matrix.median(axis=1),
        "iqr_rank": rank_matrix.quantile(0.75, axis=1) - rank_matrix.quantile(0.25, axis=1),
        "frac_runs_in_top_n": in_top.mean(axis=1),
        "worst_rank": rank_matrix.max(axis=1),
        "n_runs": rank_matrix.notna().sum(axis=1),
    }).sort_values("median_rank")


def main(sm):
    log = setup_logging(sm.log[0])
    ranking = pd.read_csv(sm.input[0], index_col=0)
    runs = sorted(__import__("glob").glob("results/club/resample/*/club_tf_ranking.csv"))
    if not runs:
        log.warning("no resampling runs found; writing the point estimate only. "
                    "Robustness is UNASSESSED - do not report stability from this file.")
        out = ranking[["rank"]].rename(columns={"rank": "median_rank"})
        out["frac_runs_in_top_n"] = np.nan
        out["n_runs"] = 0
        out.to_csv(sm.output[0])
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "robustness not assessed\n(no resampling runs)",
                ha="center", va="center")
        ax.axis("off")
        fig.savefig(sm.output[1], dpi=150, bbox_inches="tight")
        return

    mat = pd.DataFrame({r: pd.read_csv(r, index_col=0)["rank"] for r in runs})
    stability = summarize_stability(mat)
    stability.to_csv(sm.output[0])

    top = stability.head(25).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.barh(top.index, top["frac_runs_in_top_n"])
    ax.set_xlabel(f"fraction of {mat.shape[1]} resampling runs with the TF in the top 25")
    ax.set_xlim(0, 1)
    fig.savefig(sm.output[1], dpi=200, bbox_inches="tight")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
