"""Hyperparameter sweep over the gating contrasts.

Scope of the sweep is deliberately limited. Everything in the grid is downstream
of integration, so a configuration can be evaluated from the annotated atlas
without re-integrating - and, more importantly, the sweep cannot change the
composition of the atlas it is being scored on.

What the sweep may see: the gating contrasts (basal->goblet, basal->ciliated) and
nothing else. `gating_recovery_score` is computed from `score_gating` output only.
The club contrast is not computed here, and the rare-lineage contrasts are not
either. This is the mechanism behind the claim in the README that the controls
were not tuned against the answer.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import pandas as pd
import scanpy as sc

from basal_to_club.benchmark.score import gating_recovery_score, score_gating
from basal_to_club.secretory.resolver import classify, fit_thresholds, score_programs
from basal_to_club.signature.pseudobulk import (
    limma_voom_signature,
    make_pseudobulk,
    zscore_signature,
)
from basal_to_club.utils.io import load_yaml, setup_logging
from basal_to_club.viper.run_viper import fdr, load_regulon, run_pyviper


def expand_grid(grid: dict, max_configs: int) -> list[dict]:
    """Cartesian product of the grid, as a list of dotted-key -> value dicts."""
    keys = sorted(grid)
    combos = list(itertools.product(*(grid[k] for k in keys)))
    if len(combos) > max_configs:
        raise ValueError(
            f"sweep grid expands to {len(combos)} configurations, above max_configs="
            f"{max_configs}. Narrow the grid rather than raising the cap: every "
            f"configuration is a full VIPER run."
        )
    return [dict(zip(keys, values, strict=True)) for values in combos]


def apply_overrides(params: dict, overrides: dict) -> dict:
    """Return a deep-ish copy of `params` with dotted-key overrides applied."""
    import copy
    out = copy.deepcopy(params)
    for dotted, value in overrides.items():
        section, key = dotted.split(".", 1)
        out.setdefault(section, {})[key] = value
    return out


def evaluate(adata, regulons: dict[str, pd.DataFrame], params: dict,
             markers: dict, bench: dict, log) -> dict:
    """Resolve -> pseudobulk -> signature -> VIPER -> gating score, one config."""
    sec = params["secretory"]
    scores = score_programs(adata, markers)
    thresholds = fit_thresholds(scores, min_separation=sec["min_component_separation"])
    adata.obs["cell_class"] = classify(scores, thresholds,
                                       hybrid_class=sec.get("hybrid_class", True))

    sig_cfg = params["signature"]
    counts, meta = make_pseudobulk(adata, group_by=["donor_id", "cell_class"],
                                   min_cells=sig_cfg["min_cells_per_pseudobulk"])

    activity = {}
    for family, regulon in regulons.items():
        frames = []
        for contrast in bench["gating"]["contrasts"]:
            source, target = contrast.replace("basal_to_", "basal|").split("|")
            res = limma_voom_signature(counts, meta, group_col="cell_class",
                                       contrast=(target, source))
            ges = zscore_signature(res["t_stat"])
            nes = run_pyviper(ges.to_frame(contrast), regulon,
                              n_perm=params["viper"]["n_perm"],
                              min_regulon_size=params["viper"]["min_regulon_size"])
            frame = nes.rename(columns={nes.columns[0]: "nes"}).reset_index()
            frame.columns = ["tf", "nes"]
            frame["contrast"] = contrast
            frame["fdr"] = fdr(frame["nes"])
            frames.append(frame)
        activity[family] = pd.concat(frames, ignore_index=True)

    scores_df, gate = score_gating(activity, bench)
    return {"gating_recovery_score": gating_recovery_score(scores_df),
            "gate_passed": bool(gate["passed"]),
            "specificity_passed": bool(gate.get("specificity_passed", False)),
            "per_contrast": gate["per_contrast"]}


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    markers = load_yaml(sm.input.markers)["markers"]
    bench = load_yaml(sm.input.benchmark)

    cfg = params["sweep"]
    if not cfg.get("enabled", True):
        raise RuntimeError(
            "sweep.enabled is false but the sweep target was requested. Either "
            "enable it or freeze a configuration by hand and record why.")

    configs = expand_grid(cfg["grid"], cfg.get("max_configs", 24))
    log.info("sweep: %d configurations over %s", len(configs), sorted(cfg["grid"]))

    adata = sc.read_h5ad(sm.input.h5ad)
    regulons = {family: load_regulon(path)
                for family, path in zip(sm.params.families, sm.input.networks, strict=True)}

    rows = []
    for i, overrides in enumerate(configs):
        cfg_id = f"cfg{i:03d}"
        result = evaluate(adata.copy(), regulons, apply_overrides(params, overrides),
                          markers, bench, log)
        rows.append({"config_id": cfg_id, **overrides,
                     "gating_recovery_score": result["gating_recovery_score"],
                     "gate_passed": result["gate_passed"],
                     "specificity_passed": result["specificity_passed"]})
        log.info("%s score=%.4f passed=%s %s", cfg_id, result["gating_recovery_score"],
                 result["gate_passed"], overrides)
        Path(sm.output.details).parent.mkdir(parents=True, exist_ok=True)
        with open(Path(sm.output.details).parent / f"{cfg_id}.json", "w") as fh:
            json.dump({"overrides": overrides, **result}, fh, indent=2, default=str)

    out = pd.DataFrame(rows).sort_values("gating_recovery_score", ascending=False)
    out.to_csv(sm.output.table, index=False)
    with open(sm.output.details, "w") as fh:
        json.dump({"n_configs": len(rows), "grid": cfg["grid"],
                   "scored_on": sorted(bench["gating"]["contrasts"]),
                   "note": "club and rare-lineage contrasts were not computed"},
                  fh, indent=2)
    log.info("sweep complete: best %s at %.4f", out.iloc[0]["config_id"],
             out.iloc[0]["gating_recovery_score"])


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
