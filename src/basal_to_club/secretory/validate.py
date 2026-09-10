"""Two independent validations of the resolver, both against GSE121600.

1. Held-out in vivo bronchial biopsy (GSM3439925): balanced accuracy >= 0.85
   against author labels the resolver never saw.
2. ALI time course: club fraction must rise with ALI day, and the BEGM vs
   PneumaCult comparison must move the goblet fraction. A resolver that passes
   the first and fails the second is cutting a continuum at an arbitrary point.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from basal_to_club.secretory.resolver import balanced_accuracy, monotonic_club_fraction
from basal_to_club.utils.io import load_yaml, setup_logging, write_json


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    adata = sc.read_h5ad(sm.input.h5ad)
    obs = adata.obs
    holdout_id = params["secretory"]["holdout_sample"]

    mask = obs["sample_id"].astype(str) == holdout_id
    if mask.sum() == 0:
        raise ValueError(f"holdout sample {holdout_id} absent from the atlas")
    truth = obs.loc[mask, "author_cell_type"].astype(str).str.lower().replace(
        {"secretory": "club", "club cell": "club", "goblet cell": "goblet"})
    bacc = balanced_accuracy(truth, obs.loc[mask, "cell_class"].astype(str))

    ali = obs[obs["dataset_id"].astype(str) == "ruizgarcia_2019"]
    timecourse = monotonic_club_fraction(ali)

    media = {}
    if "media" in ali and ali["media"].notna().any():
        media = (ali.groupby("media", observed=True)["cell_class"]
                    .apply(lambda s: {"club": float((s == "club").mean()),
                                      "goblet": float((s == "goblet").mean()),
                                      "n": int(len(s))}).to_dict())

    threshold = params["secretory"]["holdout_min_balanced_accuracy"]
    metrics = {"holdout_sample": holdout_id, "n_holdout_cells": int(mask.sum()),
               "balanced_accuracy": bacc, "threshold": threshold,
               "holdout_passed": bool(pd.notna(bacc) and bacc >= threshold),
               "ali_timecourse": timecourse, "composition_by_media": media}
    log.info("resolver holdout balanced accuracy = %s (threshold %s)", bacc, threshold)
    if not metrics["holdout_passed"]:
        log.error("resolver failed its holdout. Downstream club results are not "
                  "interpretable until this is addressed.")
    write_json(sm.output.metrics, metrics)

    fig, ax = plt.subplots(figsize=(6, 4))
    frac = timecourse.get("club_fraction_by_day", {})
    if frac:
        days = sorted(frac)
        ax.plot(days, [frac[d] for d in days], marker="o")
    ax.set_xlabel("ALI day (GSE121600)")
    ax.set_ylabel("club fraction of secretory cells")
    ax.set_title("Club abundance over ALI differentiation")
    fig.savefig(sm.output.fig, dpi=200, bbox_inches="tight")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
