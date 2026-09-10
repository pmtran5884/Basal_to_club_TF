"""Coarse annotation: score every marker panel, arbitrate, flag disagreement."""
from __future__ import annotations

import pandas as pd
import scanpy as sc

from basal_to_club.annotate.arbitrate import marker_call, modal_of_three
from basal_to_club.utils.genes import match_panel
from basal_to_club.utils.io import load_yaml, setup_logging


def score_all_panels(adata, markers: dict) -> pd.DataFrame:
    scores, missing = {}, {}
    for name, panel in markers["markers"].items():
        if name in ("immune", "endothelial", "stromal"):
            continue
        genes, absent = match_panel(panel, adata.var_names)
        missing[name] = absent
        if not genes:
            continue
        sc.tl.score_genes(adata, genes, score_name=f"_s_{name}")
        scores[name] = adata.obs.pop(f"_s_{name}")
    adata.uns["marker_panel_missing_genes"] = missing
    return pd.DataFrame(scores, index=adata.obs_names)


def main(sm):
    log = setup_logging(sm.log[0])
    params, markers = load_yaml("config/params.yaml"), load_yaml(sm.input.markers)
    adata = sc.read_h5ad(sm.input.h5ad)

    scores = score_all_panels(adata, markers)
    # Score by cluster, not by cell: single-cell panel scores are dropout-noisy.
    by_cluster = scores.groupby(adata.obs["leiden"], observed=True).mean()
    cluster_call = marker_call(by_cluster)
    adata.obs["marker_label"] = adata.obs["leiden"].map(cluster_call).astype(str)

    consensus, disagreement = modal_of_three(
        adata.obs.get("author_cell_type", pd.Series(pd.NA, index=adata.obs_names)),
        adata.obs.get("transferred_label", pd.Series(pd.NA, index=adata.obs_names)),
        adata.obs["marker_label"],
        adata.obs.get("transfer_confidence"),
        params["annotate"]["min_transfer_confidence"])
    adata.obs["cell_type_coarse"] = consensus.astype("category")
    for name in scores.columns:
        adata.obs[f"score_{name}"] = scores[name]

    log.info("annotation disagreement on %d/%d cells (%.1f%%)", len(disagreement),
             adata.n_obs, 100 * len(disagreement) / max(adata.n_obs, 1))
    disagreement.to_csv(sm.output.disagreement)
    adata.write_h5ad(sm.output.h5ad, compression="gzip")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
