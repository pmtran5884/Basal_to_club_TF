"""Basal -> club trajectory.

Goblet, ionocyte, tuft and PNEC branches are excluded from the fate computation
so the club signature is not contaminated by cells committed elsewhere. In
GSE121600 the trajectory is anchored on real ALI day rather than pseudotime;
elsewhere pseudotime is used and flagged as such.
"""
from __future__ import annotations

import pandas as pd
import scanpy as sc

from basal_to_club.utils.io import load_yaml, setup_logging


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    cfg = params["trajectory"]
    adata = sc.read_h5ad(sm.input.h5ad)

    excluded = set(cfg.get("exclude_branches", []))
    keep = ~adata.obs["cell_class"].astype(str).isin(excluded)
    sub = adata[keep].copy()
    log.info("trajectory on %d cells (excluded branches: %s)", sub.n_obs, sorted(excluded))

    try:
        import cellrank as cr
        sc.tl.diffmap(sub)
        root = sub.obs["cell_class"].astype(str).eq("basal")
        sub.uns["iroot"] = int(root.to_numpy().nonzero()[0][0])
        sc.tl.dpt(sub)
        kernel = cr.kernels.PseudotimeKernel(sub, time_key="dpt_pseudotime").compute_transition_matrix()
        estimator = cr.estimators.GPCCA(kernel)
        estimator.compute_macrostates(n_states=6, cluster_key="cell_class")
        estimator.set_terminal_states_from_macrostates()
        estimator.compute_fate_probabilities()
        fates = pd.DataFrame(estimator.fate_probabilities, index=sub.obs_names,
                             columns=[str(c) for c in estimator.terminal_states.cat.categories])
    except Exception as exc:                            # noqa: BLE001
        log.warning("CellRank unavailable (%s); falling back to DPT pseudotime only", exc)
        fates = pd.DataFrame({"dpt_pseudotime": sub.obs.get("dpt_pseudotime", pd.NA)},
                             index=sub.obs_names)

    # Where real differentiation time exists, use it: GSE121600 ALI day.
    if "ali_day" in sub.obs and sub.obs["ali_day"].notna().any():
        fates["ali_day"] = sub.obs["ali_day"]
        fates["time_anchor"] = "observed_ali_day"
    else:
        fates["time_anchor"] = "pseudotime"

    fates.to_csv(sm.output.fates)
    sub.obsm["fate_probabilities"] = fates.select_dtypes("number").to_numpy()
    sub.write_h5ad(sm.output.h5ad, compression="gzip")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
