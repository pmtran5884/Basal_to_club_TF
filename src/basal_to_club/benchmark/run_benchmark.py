"""Score the frozen positive controls and write the gate decision."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from basal_to_club.benchmark.score import score_gating, score_reported_only
from basal_to_club.utils.io import load_yaml, setup_logging, write_json


def lineage_cell_counts(adata) -> dict[str, dict]:
    obs = adata.obs
    out = {}
    for lineage, grp in obs.groupby(obs["cell_class"].astype(str), observed=True):
        out[lineage] = {"n_cells": int(len(grp)),
                        "n_donors": int(grp["donor_id"].nunique())}
    return out


def main(sm):
    log = setup_logging(sm.log[0])
    bench = load_yaml(sm.input.bench)
    by_family = pd.read_csv(sm.input.by_family)
    combined = pd.read_csv(sm.input.activity)
    counts = lineage_cell_counts(sc.read_h5ad(sm.input.atlas))

    families = {fam: sub for fam, sub in by_family.groupby("family", observed=True)}
    scores, gate = score_gating(families, bench)
    rare = score_reported_only(combined, bench, counts)

    gate["lineage_cell_counts"] = counts
    gate["reported_only"] = rare.to_dict("records")
    gate["frozen_config"] = bool(load_yaml("config/params.yaml").get("_frozen", False))

    pd.concat([scores, rare.assign(kind="reported_only")], ignore_index=True) \
      .to_csv(sm.output.scores, index=False)
    write_json(sm.output.gate, gate)
    log.info("GATE %s | recovery score %.3f",
             "PASSED" if gate["passed"] else "FAILED", gate["gating_recovery_score"])
    for contrast, detail in gate["per_contrast"].items():
        log.info("  %-24s %s %s", contrast, "pass" if detail["passed"] else "FAIL",
                 detail["failure_reasons"] or "")

    prim = scores[(scores["kind"] == "primary") & scores["rank_percentile"].notna()]
    fig, ax = plt.subplots(figsize=(7, 4))
    for fam, sub in prim.groupby("family", observed=True):
        ax.scatter(sub["contrast"], 100 * (1 - sub["rank_percentile"]), label=fam, s=70)
    if len(rare):
        ev = rare[rare["rank_percentile"].notna()]
        if len(ev):
            ax.scatter(ev["contrast"], 100 * (1 - ev["rank_percentile"]),
                       marker="^", facecolors="none", edgecolors="grey",
                       label="rare (reported, not gating)")
    ax.axhline(95, ls="--", lw=1, color="k")
    ax.set_ylabel("recovery percentile of the lineage TF")
    ax.set_ylim(0, 101)
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(sm.output.curve, dpi=200, bbox_inches="tight")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
