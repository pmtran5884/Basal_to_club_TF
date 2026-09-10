"""Per-sample adaptive QC.

Thresholds are MAD-based within sample rather than global: a global mito cut
preferentially deletes ciliated cells, which are genuinely mitochondria-rich, and
that would bias the very lineage the benchmark depends on.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from basal_to_club.utils.io import load_yaml, setup_logging, write_json


def mad_outlier(values: np.ndarray, n_mads: float = 3.0, side: str = "both") -> np.ndarray:
    med = np.median(values)
    mad = np.median(np.abs(values - med)) or 1e-9
    lo, hi = med - n_mads * 1.4826 * mad, med + n_mads * 1.4826 * mad
    if side == "upper":
        return values > hi
    if side == "lower":
        return values < lo
    return (values < lo) | (values > hi)


def qc_metrics(adata) -> None:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)


def flag_low_quality(adata, n_mads: float, min_genes: int, max_pct_mito: float,
                     sample_key: str = "sample_id") -> pd.Series:
    flags = pd.Series(False, index=adata.obs_names)
    for _, idx in adata.obs.groupby(sample_key, observed=True).groups.items():
        o = adata.obs.loc[idx]
        bad = (mad_outlier(o["log1p_total_counts"].to_numpy(), n_mads)
               | mad_outlier(o["log1p_n_genes_by_counts"].to_numpy(), n_mads)
               | mad_outlier(o["pct_counts_mt"].to_numpy(), n_mads, side="upper")
               | (o["n_genes_by_counts"].to_numpy() < min_genes)
               | (o["pct_counts_mt"].to_numpy() > max_pct_mito))
        flags.loc[idx] = bad
    return flags


def drop_nonepithelial(adata, markers: dict, min_score: float = 0.1):
    """Remove immune/endothelial/stromal cells. Epithelial-only from here on."""
    from basal_to_club.utils.genes import match_panel
    keep = pd.Series(True, index=adata.obs_names)
    for lineage in ("immune", "endothelial", "stromal"):
        genes, _ = match_panel(markers["markers"][lineage], adata.var_names)
        if not genes:
            continue
        sc.tl.score_genes(adata, genes, score_name="_nonep")
        keep &= adata.obs.pop("_nonep") < min_score
    return adata[keep].copy()


def main(sm):
    log = setup_logging(sm.log[0])
    params, markers = load_yaml("config/params.yaml"), load_yaml("config/markers.yaml")
    q = params["qc"]
    adata = sc.read_h5ad(sm.input.h5ad)
    n0 = adata.n_obs
    qc_metrics(adata)

    if q.get("doublets") == "scdblfinder":
        try:
            import scanpy.external as sce
            sce.pp.scrublet(adata, batch_key="sample_id")
            adata = adata[~adata.obs["predicted_doublet"]].copy()
        except Exception as exc:                       # noqa: BLE001
            log.warning("doublet detection unavailable (%s); continuing and recording it", exc)

    adata = adata[~flag_low_quality(adata, q["n_mads"], q["min_genes"], q["max_pct_mito"])].copy()
    if q.get("drop_nonepithelial", True):
        adata = drop_nonepithelial(adata, markers)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.write_h5ad(sm.output.h5ad, compression="gzip")

    write_json(sm.output.metrics, {
        "dataset_id": str(adata.obs["dataset_id"].iloc[0]), "n_cells_in": int(n0),
        "n_cells_out": int(adata.n_obs), "frac_retained": float(adata.n_obs / max(n0, 1)),
        "median_genes": float(adata.obs["n_genes_by_counts"].median()),
        "median_pct_mt": float(adata.obs["pct_counts_mt"].median()),
        "n_samples": int(adata.obs["sample_id"].nunique())})
    fig = sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                       groupby="sample_id", rotation=90, show=False, return_fig=True)
    fig.savefig(sm.output.fig, dpi=150, bbox_inches="tight")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
