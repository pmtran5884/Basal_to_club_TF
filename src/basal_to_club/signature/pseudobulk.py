"""Pseudobulk construction and differential signatures.

Donor is the unit of replication, not the cell. Cell-level differential tests
treat thousands of cells from one donor as independent and produce p-values that
are essentially a function of cell count; a VIPER run on that signature inherits
the inflation directly.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def make_pseudobulk(adata, group_by: list[str], min_cells: int = 25,
                    layer: str = "counts") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum counts within (dataset, donor, cell_class). Returns (counts, metadata)."""
    import scanpy as sc

    obs = adata.obs
    key = obs[group_by].astype(str).agg("|".join, axis=1)
    sizes = key.value_counts()
    keep_groups = sizes[sizes >= min_cells].index
    mask = key.isin(keep_groups)
    if mask.sum() == 0:
        raise ValueError(f"no pseudobulk group reached min_cells={min_cells}")

    # scanpy's aggregate takes the NAME of an obs column, not an array of labels;
    # passing the array raises KeyError on the first label. Stage it as a
    # temporary column and remove it again so the caller's object is unchanged.
    tmp = "_pb_group_key"
    adata.obs[tmp] = pd.Categorical(key)
    try:
        agg = sc.get.aggregate(adata[mask], by=tmp, func="sum", layer=layer)
    finally:
        del adata.obs[tmp]
    counts = pd.DataFrame(np.asarray(agg.layers["sum"]), index=agg.obs_names,
                          columns=agg.var_names)
    # strict=True: a group key that does not split into exactly len(group_by)
    # fields means a "|" in a donor or cell-class name has corrupted the grouping.
    meta = pd.DataFrame([dict(zip(group_by, g.split("|"), strict=True)) for g in counts.index],
                        index=counts.index)
    meta["n_cells"] = sizes.reindex(counts.index).to_numpy()
    return counts, meta


def limma_voom_signature(counts: pd.DataFrame, meta: pd.DataFrame, group_col: str,
                         case: str, control: str, block: str = "donor_id",
                         covariates: list[str] | None = None) -> pd.DataFrame:
    """limma-voom with donor blocking. Returns a t-statistic GES per gene.

    VIPER consumes the t-statistic, not the log fold change: the t-statistic
    already carries the within-group variance, which is what a rank-based
    enrichment method needs.
    """
    from rpy2.robjects import default_converter, pandas2ri, r
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    # pandas2ri.activate() was removed in rpy2 3.6; conversion is scoped to a
    # context manager now. The global-activation form raises DeprecationWarning
    # as an exception, so it fails loudly rather than degrading.
    cv = default_converter + pandas2ri.converter
    # Loaded for the side effect of attaching the packages in the R namespace;
    # the model is fitted through r() strings below, not through these handles.
    importr("limma"), importr("edgeR")
    sub = meta[meta[group_col].isin([case, control])].copy()
    y = counts.loc[sub.index].T

    # Direction is set explicitly from `case`/`control`. Deriving it from column
    # order instead would silently invert every contrast whose control sorts
    # first alphabetically (e.g. basal vs club), putting the lineage TF at the
    # BOTTOM of the ranking; see tests/test_signature_direction.py.
    sub["grp_contrast"] = np.where(sub[group_col].to_numpy() == case, "case", "control")

    with localconverter(cv):
        r.assign("counts", y)
        r.assign("meta", sub.reset_index(drop=True))
    terms = " + ".join([f"factor({c})" for c in (covariates or []) if sub[c].nunique() > 1])
    formula = "~0 + factor(grp_contrast)" + (f" + {terms}" if terms else "")
    r(f'''
    dge <- edgeR::DGEList(counts = as.matrix(counts))
    keep <- edgeR::filterByExpr(dge, group = meta$grp_contrast)
    dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
    design <- model.matrix({formula}, data = meta)
    colnames(design) <- make.names(colnames(design))
    v <- limma::voom(dge, design)
    corfit <- limma::duplicateCorrelation(v, design, block = meta${block})
    v <- limma::voom(dge, design, block = meta${block}, correlation = corfit$consensus)
    fit <- limma::lmFit(v, design, block = meta${block}, correlation = corfit$consensus)
    cm <- limma::makeContrasts(contrasts = "factor.grp_contrast.case-factor.grp_contrast.control",
                               levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
    res <- limma::topTable(fit2, number = Inf, sort.by = "none")
    ''')
    with localconverter(cv):
        res = pd.DataFrame(r("res"))
    res.index = list(r("rownames(res)"))
    return res.rename(columns={"t": "t_stat", "logFC": "log_fc", "adj.P.Val": "fdr"})


def zscore_signature(t_stats: pd.Series) -> pd.Series:
    """VIPER expects a z-like GES. Rank-normalize, then map to the normal quantile."""
    from scipy.stats import norm
    ranks = t_stats.rank(method="average") / (len(t_stats) + 1)
    return pd.Series(norm.ppf(ranks), index=t_stats.index, name="ges")
