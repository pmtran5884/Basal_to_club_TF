"""metaVIPER integration across networks, and per-family aggregation.

Families are kept separate as well as combined, because the benchmark requires a
control TF to pass in BOTH the de novo and the prior family. Agreement between an
expression-inferred and a prior-knowledge network is the evidence that a result is
not an artifact of one network construction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from basal_to_club.utils.io import setup_logging


def family_of(network: str) -> str:
    return "de_novo" if network.startswith("aracne") else "prior"


def integrate(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Stouffer combination across networks, weighted by |NES| as metaVIPER does.

    A network with no opinion about a TF contributes little rather than pulling
    the estimate toward zero.
    """
    def combine(g):
        w = np.abs(g["nes"].to_numpy())
        w = w / w.sum() if w.sum() > 0 else np.full(len(g), 1 / len(g))
        return pd.Series({"nes": float((g["nes"].to_numpy() * w).sum() * np.sqrt(len(g))),
                          "n_networks": int(len(g)),
                          "fdr": float(g["fdr"].min())})
    return df.groupby(group_cols, observed=True).apply(combine, include_groups=False).reset_index()


def main(sm):
    log = setup_logging(sm.log[0])
    parts = [pd.read_csv(p) for p in sm.input]
    allact = pd.concat(parts, ignore_index=True)
    allact["family"] = allact["network"].map(family_of)

    combined = integrate(allact, ["tf", "contrast", "method"])
    by_family = integrate(allact, ["tf", "contrast", "method", "family"])
    log.info("metaVIPER: %d TFs across %d contrasts, %d networks",
             combined["tf"].nunique(), combined["contrast"].nunique(),
             allact["network"].nunique())
    combined.to_csv(sm.output.combined, index=False)
    by_family.to_csv(sm.output.by_family, index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
