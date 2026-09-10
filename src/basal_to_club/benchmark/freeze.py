"""Freeze the winning sweep configuration.

The sweep is scored on `gating_recovery_score`, which by construction cannot see
the rare-lineage controls or the club contrast. Once frozen, config/frozen.yaml
overrides params.yaml for every subsequent run, so the club prediction is made
under a configuration that was selected before anyone looked at a club result.
"""
from __future__ import annotations

import pandas as pd
import yaml

from basal_to_club.utils.io import setup_logging


def main(sm):
    log = setup_logging(sm.log[0])
    sweep = pd.read_csv(sm.input[0])
    best = sweep.sort_values("gating_recovery_score", ascending=False).iloc[0]
    frozen = {k.split("param.")[1]: yaml.safe_load(str(v))
              for k, v in best.items() if k.startswith("param.")}

    nested: dict = {}
    for dotted, value in frozen.items():
        node = nested
        *parents, leaf = dotted.split(".")
        for p in parents:
            node = node.setdefault(p, {})
        node[leaf] = value
    nested["_frozen_from"] = {"gating_recovery_score": float(best["gating_recovery_score"]),
                              "sweep_row": int(best.name)}
    with open(sm.output[0], "w") as fh:
        yaml.safe_dump(nested, fh, sort_keys=True)
    log.info("froze configuration with gating recovery score %.4f",
             best["gating_recovery_score"])


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
