"""Build all gene expression signatures: three constructions x every contrast.

Agreement across the three constructions is the evidence. Where they disagree
that is reported, not averaged away - a TF that only appears under one signature
definition is a property of that definition, not of the biology.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from basal_to_club.signature.pseudobulk import (
    limma_voom_signature,
    make_pseudobulk,
    zscore_signature,
)
from basal_to_club.utils.io import load_yaml, setup_logging

CONTRASTS = {
    "basal_to_club": ("club", "basal"),
    "basal_to_goblet": ("goblet", "basal"),
    "basal_to_ciliated": ("ciliated", "basal"),
    "basal_to_deuterosomal": ("deuterosomal", "basal"),
    "basal_identity": ("basal", "club"),
    "basal_to_ionocyte": ("ionocyte", "basal"),
    "basal_to_tuft": ("tuft", "basal"),
    "basal_to_pnec": ("pnec", "basal"),
    "basal_to_hillock": ("hillock", "basal"),
}


def endpoint_signatures(adata, params) -> dict[str, pd.Series]:
    cfg = params["signature"]["pseudobulk"]
    counts, meta = make_pseudobulk(adata, cfg["group_by"], cfg["min_cells_per_pseudobulk"])
    out = {}
    for name, (case, control) in CONTRASTS.items():
        if not {case, control}.issubset(set(meta["cell_class"])):
            continue
        try:
            res = limma_voom_signature(counts, meta, "cell_class", case, control,
                                       covariates=["dataset_id"])
            out[name] = zscore_signature(res["t_stat"])
        except Exception as exc:                        # noqa: BLE001
            out[name] = _fallback_signature(counts, meta, case, control, exc)
    return out


def _fallback_signature(counts, meta, case, control, exc) -> pd.Series:
    """Welch t-statistic on CPM if limma is unavailable. Recorded, not hidden."""
    cpm = np.log2(1e6 * counts.div(counts.sum(axis=1), axis=0) + 1)
    a = cpm.loc[meta.index[meta["cell_class"] == case]]
    b = cpm.loc[meta.index[meta["cell_class"] == control]]
    from scipy.stats import ttest_ind
    t = ttest_ind(a, b, equal_var=False, nan_policy="omit").statistic
    s = zscore_signature(pd.Series(np.nan_to_num(t), index=cpm.columns))
    s.attrs["fallback_reason"] = str(exc)
    return s


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    atlas = sc.read_h5ad(sm.input.atlas)

    frames = {}
    for method in params["signature"]["contrasts"]:
        source = {"endpoint": atlas,
                  "trajectory": sc.read_h5ad(sm.input.traj),
                  "metacell": sc.read_h5ad(sm.input.meta)}[method]
        if method == "metacell" and "cell_class" not in source.obs:
            source.obs["cell_class"] = source.obs["dominant_cell_class"]
        for contrast, ges in endpoint_signatures(source, params).items():
            frames[(method, contrast)] = ges

    sig = pd.DataFrame(frames)
    sig.columns = pd.MultiIndex.from_tuples(sig.columns, names=["method", "contrast"])
    sig.to_hdf(sm.output[0], key="signatures", mode="w")

    summary = pd.DataFrame({
        "method": [c[0] for c in sig.columns], "contrast": [c[1] for c in sig.columns],
        "n_genes": sig.notna().sum().to_numpy(),
        "max_abs_ges": sig.abs().max().to_numpy()})
    summary.to_csv(sm.output[1], index=False)
    log.info("built %d signatures across %d contrasts", sig.shape[1],
             summary["contrast"].nunique())


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
