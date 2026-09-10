"""Turn one ChEA3 library into a signed, weighted regulon."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from basal_to_club.network.prior import (
    audit_signs_against_reference,
    build_regulon,
    cross_library_support,
    read_gmt,
)
from basal_to_club.utils.io import load_yaml, setup_logging, write_json


def main(sm):
    log = setup_logging(sm.log[0])
    params, regcfg = load_yaml("config/params.yaml"), load_yaml("config/regulons.yaml")
    chea = regcfg["prior"]["chea3"]
    libdir = Path(sm.input.libs)

    libraries = {lib: read_gmt(libdir / f"{lib}.gmt") for lib in chea["libraries"]}
    support = cross_library_support(libraries)

    mc = sc.read_h5ad(sm.input.meta)
    expr = pd.DataFrame(np.asarray(mc.X.todense() if hasattr(mc.X, "todense") else mc.X),
                        index=mc.obs_names, columns=[str(v).upper() for v in mc.var_names])

    regulon = build_regulon(libraries[sm.wildcards.library], expr, support,
                            min_abs_rho=chea["min_abs_rho"],
                            min_targets=params["regulon_filter"]["min_targets"],
                            max_targets=params["regulon_filter"]["max_targets"])
    log.info("chea3 %s -> %d TFs, %d edges, %.1f%% edges with likelihood 0",
             sm.wildcards.library, regulon["tf"].nunique(), len(regulon),
             100 * float((regulon["likelihood"] == 0).mean()))

    ref_path = Path("results/networks/collectri.tsv")
    if ref_path.exists():
        audit = audit_signs_against_reference(regulon, pd.read_csv(ref_path, sep="\t"))
        write_json(f"results/networks/{sm.wildcards.library}_sign_audit.json", audit)
        log.info("sign audit vs CollecTRI: %s", audit)
    regulon.to_csv(sm.output[0], sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
