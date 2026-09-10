"""Gate every network before VIPER sees it.

A network where FOXJ1 does not regulate cilia genes is broken in a way that no
downstream statistic makes obvious, so it is caught here rather than explained
away later.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from basal_to_club.network.prior import regulon_sanity_check
from basal_to_club.utils.io import load_yaml, setup_logging


def main(sm):
    log = setup_logging(sm.log[0])
    rules = load_yaml("config/regulons.yaml")["network_qc"]["required_regulon_sanity"]
    rows = []
    for path in sm.input:
        name = Path(path).stem
        reg = pd.read_csv(path, sep="\t")
        check = regulon_sanity_check(reg, rules).assign(network=name)
        rows.append(check)
        summary = {"network": name, "n_tfs": reg["tf"].nunique(), "n_edges": len(reg),
                   "median_regulon_size": float(reg.groupby("tf").size().median()),
                   "frac_repressive": float((reg["mo"] < 0).mean())}
        log.info("%s", summary)
        if not check["passed"].all():
            log.error("%s failed regulon sanity: %s", name,
                      check.loc[~check["passed"], "tf"].tolist())
    pd.concat(rows, ignore_index=True).to_csv(sm.output[0], index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
