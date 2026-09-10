"""Raw files -> AnnData with a harmonized obs schema.

Sample-level exclusions are asserted HERE and fail loudly. GSE121600 contains
mouse and pig samples; GSE166766 contains SARS-CoV-2-infected cells whose
interferon program would dominate any differential signature. Leaving either to
a downstream subset is how they end up in the analysis.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from basal_to_club.utils.io import setup_logging, write_provenance

OBS_SCHEMA = ["dataset_id", "sample_id", "donor_id", "system", "tissue", "chemistry",
              "media", "ali_day", "timepoint", "condition", "airway_position",
              "smoking_status", "author_cell_type"]


def harmonize_obs(adata: ad.AnnData, spec: dict) -> ad.AnnData:
    """Impose the common obs schema; missing fields become NA rather than absent."""
    obs = adata.obs
    obs["dataset_id"] = spec["id"]
    for col in OBS_SCHEMA:
        if col not in obs:
            obs[col] = pd.NA
    per_sample = spec.get("samples") or {}
    if per_sample and "sample_id" in obs:
        meta = pd.DataFrame.from_dict(per_sample, orient="index")
        for col in meta.columns:
            mapped = obs["sample_id"].map(meta[col])
            obs[col] = obs[col].where(mapped.isna(), mapped)
    adata.obs = obs[OBS_SCHEMA + [c for c in obs.columns if c not in OBS_SCHEMA]]
    return adata


def apply_exclusions(adata: ad.AnnData, spec: dict, log) -> ad.AnnData:
    n0 = adata.n_obs
    drop = set(spec.get("exclude_samples") or [])
    if drop:
        keep = ~adata.obs["sample_id"].astype(str).isin(drop)
        if keep.sum() == adata.n_obs:
            raise ValueError(
                f"{spec['id']}: exclude_samples {sorted(drop)} matched nothing. "
                "The sample id convention changed - fix the mapping rather than proceeding.")
        adata = adata[keep].copy()
        log.info("%s: dropped %d cells from excluded samples %s",
                 spec["id"], n0 - adata.n_obs, sorted(drop))

    dfilter = spec.get("disease_filter")
    if dfilter and "condition" in adata.obs:
        cond = adata.obs["condition"].astype(str).str.lower()
        keep = cond.isin([d.lower() for d in dfilter])
        if keep.sum() == 0:
            raise ValueError(
                f"{spec['id']}: disease_filter {dfilter} kept 0 cells. Inspect the "
                f"condition vocabulary: {sorted(cond.unique())[:15]}")
        adata = adata[keep].copy()
        log.info("%s: disease_filter %s kept %d/%d cells", spec["id"], dfilter, adata.n_obs, n0)

    if "organism" in adata.obs:
        keep = adata.obs["organism"].astype(str).str.contains("Homo sapiens|human", case=False)
        adata = adata[keep].copy()
    return adata


def main(sm):
    log = setup_logging(sm.log[0])
    spec = sm.params.spec
    raw = Path(f"data/raw/{spec['id']}")
    h5ads = sorted(raw.glob("*.h5ad"))
    if h5ads:
        adata = ad.read_h5ad(h5ads[0])
    else:
        import scanpy as sc
        adata = sc.read_10x_mtx(raw, var_names="gene_symbols", cache=True)
    adata = apply_exclusions(harmonize_obs(adata, spec), spec, log)
    adata.var_names_make_unique()
    adata.uns["dataset_spec"] = {k: v for k, v in spec.items() if not isinstance(v, dict)}
    adata.write_h5ad(sm.output.h5ad, compression="gzip")
    write_provenance(sm.output.prov, {
        "dataset_id": spec["id"], "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "samples": sorted(adata.obs["sample_id"].astype(str).unique().tolist())})


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
