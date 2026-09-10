"""Run VIPER for one network across every signature.

Both a Python (pyviper) and an R (viper) implementation are supported; the
R package is the reference implementation and pyviper is used for throughput.
Which one produced a result is recorded in the output, because "we ran VIPER"
is not a reproducible statement on its own.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from basal_to_club.utils.io import load_yaml, setup_logging


def load_regulon(path: str) -> pd.DataFrame:
    reg = pd.read_csv(path, sep="\t")
    missing = {"tf", "target", "mo", "likelihood"} - set(reg.columns)
    if missing:
        raise ValueError(f"{path}: regulon missing columns {sorted(missing)}")
    return reg


def apply_per_signature(signatures: pd.DataFrame,
                        score_fn: "callable") -> pd.DataFrame:
    """Score each signature column separately, on its own measured genes.

    Contrasts do not share a gene universe: filterByExpr runs per contrast, so a
    gene expressed in the goblet comparison may be absent from the ciliated one.
    Merging the columns therefore produces a matrix with NaNs, and aREA
    propagates a single NaN across the entire sample -- the failure mode is an
    all-NaN activity matrix returned WITHOUT an error, which looks like a result
    until someone checks. Scoring column by column keeps every contrast's full
    gene set and makes an empty signature an explicit failure.
    """
    out = []
    for col in signatures.columns:
        s = signatures[col].dropna()
        if s.empty:
            raise ValueError(f"signature '{col}' has no non-missing genes")
        act = score_fn(s.to_frame(name=col))
        out.append(act)
    combined = pd.concat(out, axis=1)
    if combined.isna().all().any():
        bad = combined.columns[combined.isna().all()].tolist()
        raise ValueError(f"VIPER returned no scores for signatures {bad}")
    return combined


def run_pyviper(signatures: pd.DataFrame, regulon: pd.DataFrame, n_perm: int,
                pleiotropy: bool, seed: int = 0) -> pd.DataFrame:
    import anndata as ad
    import pyviper

    def score_one(sig: pd.DataFrame) -> pd.DataFrame:
        interactome = pyviper.Interactome("net", regulon.rename(columns={"tf": "regulator"}))
        interactome.filter_targets(sig.index)
        ges = ad.AnnData(sig.T.to_numpy().astype(np.float32),
                         obs=pd.DataFrame(index=sig.columns),
                         var=pd.DataFrame(index=sig.index))
        act = pyviper.viper(gex_data=ges, interactome=interactome, enrichment="area",
                            njobs=1, verbose=False)
        return pd.DataFrame(act.X, index=act.obs_names, columns=act.var_names).T

    return apply_per_signature(signatures, score_one)


def run_r_viper(signatures: pd.DataFrame, regulon: pd.DataFrame, n_perm: int) -> pd.DataFrame:
    from rpy2.robjects import default_converter, pandas2ri, r
    from rpy2.robjects.conversion import localconverter
    # activate() was removed in rpy2 3.6 and now raises; conversion is scoped.
    cv = default_converter + pandas2ri.converter
    with localconverter(cv):
        r.assign("ges", signatures)
        r.assign("net", regulon)
    r('''
    suppressMessages(library(viper))
    regul <- split(net, net$tf)
    regul <- lapply(regul, function(d) {
      tfmode <- setNames(as.numeric(d$mo), d$target)
      list(tfmode = tfmode, likelihood = as.numeric(d$likelihood))
    })
    class(regul) <- "regulon"
    res <- viper(as.matrix(ges), regul, method = "none", minsize = 20, verbose = FALSE)
    ''')
    with localconverter(cv):
        return pd.DataFrame(r("as.data.frame(res)"))


def fdr(nes: pd.Series) -> pd.Series:
    from scipy.stats import norm
    from statsmodels.stats.multitest import multipletests
    p = 2 * norm.sf(np.abs(nes.to_numpy()))
    return pd.Series(multipletests(p, method="fdr_bh")[1], index=nes.index)


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    cfg = params["viper"]
    signatures = pd.read_hdf(sm.input.signatures, key="signatures")
    regulon = load_regulon(sm.input.network)

    impl_used, activity = None, None
    for impl in cfg["implementations"]:
        try:
            activity = (run_pyviper(signatures, regulon, cfg["n_permutations"],
                                    cfg["pleiotropy_correction"], params.get("seed", 0))
                        if impl == "pyviper" else
                        run_r_viper(signatures, regulon, cfg["n_permutations"]))
            impl_used = impl
            break
        except Exception as exc:                        # noqa: BLE001
            log.warning("%s failed (%s)", impl, exc)
    if activity is None:
        raise RuntimeError("no VIPER implementation succeeded")

    long = (activity.reset_index(names="tf")
                    .melt(id_vars="tf", var_name="signature", value_name="nes"))
    long[["method", "contrast"]] = pd.DataFrame(
        [tuple(c) if isinstance(c, tuple) else (None, c) for c in long["signature"]],
        index=long.index)
    long["fdr"] = long.groupby("signature")["nes"].transform(fdr)
    long["network"] = sm.wildcards.network
    long["implementation"] = impl_used
    log.info("%s: %d TFs x %d signatures via %s", sm.wildcards.network,
             activity.shape[0], activity.shape[1], impl_used)
    long.to_csv(sm.output[0], index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
