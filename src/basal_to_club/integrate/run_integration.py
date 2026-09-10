"""Batch integration.

The corrected latent space is used for neighbours, clustering and label transfer
ONLY. Every signature downstream is recomputed from raw counts within donor -
corrected expression never reaches VIPER, because batch correction moves genes in
ways that a regulon-based method will read as regulatory signal.
"""
from __future__ import annotations

import scanpy as sc

from basal_to_club.utils.io import load_yaml, setup_logging, write_json


def integrate_scvi(adata, cfg: dict, seed: int = 0):
    import scvi
    scvi.settings.seed = seed
    adata.layers.setdefault("counts", adata.X.copy())
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3",
                                batch_key=cfg["batch_key"], layer="counts", subset=False)
    sub = adata[:, adata.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(sub, layer="counts", batch_key=cfg["batch_key"],
                                  categorical_covariate_keys=cfg.get("categorical_covariates"))
    model = scvi.model.SCVI(sub, n_latent=cfg["n_latent"], n_layers=cfg["n_layers"])
    model.train(max_epochs=cfg["max_epochs"])
    adata.obsm["X_latent"] = model.get_latent_representation()
    return adata, model


def integrate_harmony(adata, cfg: dict):
    import scanpy.external as sce
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, batch_key=cfg["batch_key"])
    sc.pp.pca(adata, n_comps=50, mask_var="highly_variable")
    sce.pp.harmony_integrate(adata, key=cfg["batch_key"])
    adata.obsm["X_latent"] = adata.obsm["X_pca_harmony"]
    return adata, None


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    cfg = params["integrate"]
    adata = sc.read_h5ad(sm.input[0])
    adata, model = (integrate_scvi(adata, cfg, params.get("seed", 0))
                    if cfg["method"] == "scvi" else integrate_harmony(adata, cfg))
    sc.pp.neighbors(adata, use_rep="X_latent")
    sc.tl.leiden(adata, resolution=cfg["leiden_resolution"], key_added="leiden")
    sc.tl.umap(adata)
    adata.write_h5ad(sm.output.h5ad, compression="gzip")

    metrics = {"method": cfg["method"], "n_cells": int(adata.n_obs),
               "n_clusters": int(adata.obs["leiden"].nunique())}
    try:
        from scib_metrics.benchmark import Benchmarker
        bm = Benchmarker(adata, batch_key=cfg["batch_key"], label_key="leiden",
                         embedding_obsm_keys=["X_latent"])
        bm.benchmark()
        metrics["scib"] = bm.get_results(min_max_scale=False).to_dict()
    except Exception as exc:                            # noqa: BLE001
        log.warning("scib metrics unavailable (%s)", exc)
        metrics["scib"] = None
    write_json(sm.output.metrics, metrics)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
