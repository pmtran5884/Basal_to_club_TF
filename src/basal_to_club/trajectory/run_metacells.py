"""Metacells for ARACNe and for rare-lineage signatures.

ARACNe needs sample-level expression variation with usable dynamic range;
single cells give it dropout instead. Metacells are built WITHIN donor so an
aggregate never mixes donors, which would manufacture correlation structure
that ARACNe would read as regulatory.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from basal_to_club.utils.io import load_yaml, setup_logging


def build_metacells(adata, cells_per_metacell: int = 75, donor_key: str = "donor_id",
                    seed: int = 0):
    import anndata as ad
    try:
        import SEACells
        have_seacells = True
    except ImportError:
        have_seacells = False

    parts = []
    for donor, idx in adata.obs.groupby(donor_key, observed=True).groups.items():
        sub = adata[idx].copy()
        n_mc = max(2, sub.n_obs // cells_per_metacell)
        if have_seacells and sub.n_obs >= 5 * cells_per_metacell:
            model = SEACells.core.SEACells(sub, build_kernel_on="X_latent",
                                           n_SEACells=n_mc, convergence_epsilon=1e-5)
            model.construct_kernel_matrix()
            model.initialize_archetypes()
            model.fit(min_iter=10, max_iter=100)
            labels = sub.obs["SEACell"].astype(str)
        else:
            # Deterministic fallback: k-means on the integrated latent space.
            from sklearn.cluster import KMeans
            km = KMeans(n_clusters=n_mc, random_state=seed, n_init=10)
            labels = ("km" + pd.Series(km.fit_predict(sub.obsm["X_latent"]),
                                       index=sub.obs_names).astype(str))
        agg = sc.get.aggregate(sub, by=labels.values, func="sum", layer="counts")
        agg.obs[donor_key] = donor
        for key in ("dataset_id", "system", "media"):
            if key in sub.obs:
                agg.obs[key] = sub.obs[key].iloc[0]
        agg.obs["dominant_cell_class"] = (
            sub.obs.groupby(labels.values, observed=True)["cell_class"]
               .agg(lambda s: s.value_counts().index[0]).reindex(agg.obs_names).values)
        agg.obs["purity"] = (
            sub.obs.groupby(labels.values, observed=True)["cell_class"]
               .agg(lambda s: float(s.value_counts(normalize=True).iloc[0]))
               .reindex(agg.obs_names).values)
        parts.append(agg)
    return ad.concat(parts, index_unique="-")


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    cfg = params["trajectory"]["metacells"]
    adata = sc.read_h5ad(sm.input.h5ad)
    mc = build_metacells(adata, cfg["cells_per_metacell"], seed=params.get("seed", 0))
    mc.layers["counts"] = mc.X.copy()
    sc.pp.normalize_total(mc, target_sum=1e6)
    sc.pp.log1p(mc)
    log.info("built %d metacells; median purity %.2f", mc.n_obs, float(np.median(mc.obs["purity"])))
    mc.write_h5ad(sm.output.h5ad, compression="gzip")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
