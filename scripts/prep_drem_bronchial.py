"""Stage 1 of the bronchial DREM arm: build the GSE233145 object and both label arms.

Data: primary human bronchial epithelial cells at air-liquid interface, days
0/3/5/7/14/21/28, two pulmonary-healthy donors (the two end-stage COPD donors in
the deposit are loaded only with --health-state COPD). This replaces the nasal
GSE121600 time course used by scripts/prep_drem_ali.py; everything downstream is
unchanged, so the two tissues are directly comparable.

Arms are the same two class definitions as the VIPER arms and the nasal DREM arm:
  atlas_label : HLCA core labels transferred by the shared classifier, the
                atlas's own Club/Goblet split kept.
  resolver    : basal and ciliated as above; the secretory compartment re-split
                by the club/goblet resolver.

The deposit's own annotation has a single undivided `Secretory` class and no
barcodes in its metadata file, so it cannot seed either arm -- it is used as an
independent composition check instead (published_vs_called_composition.csv).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from basal_to_club.data import gse233145 as G            # noqa: E402
from basal_to_club.labels import atlas_transfer as A     # noqa: E402
from basal_to_club.secretory import resolver as R        # noqa: E402

RAW = ROOT / "data/raw/gse233145"
HLCA = ROOT / "data/raw/hlca_core/hlca_core.h5ad"
OUT = ROOT / "results/drem_bronchial"
OUT.mkdir(parents=True, exist_ok=True)

MARKERS = yaml.safe_load(open(ROOT / "config/markers.yaml"))["markers"]
PARAMS = yaml.safe_load(open(ROOT / "config/params.yaml"))
N_HVG = 3000

# published class -> the coarse class it should correspond to, for the
# composition check only. Transitional Ciliated is deliberately folded into
# ciliated (it is the ascending branch), Suprabasal into basal (the atlas
# transfers cultured suprabasal/hillock cells to ann_level_3 'Basal').
PUBLISHED_TO_COARSE = {"Basal_1": "basal", "Basal_2": "basal", "Suprabasal": "basal",
                       "Secretory": "secretory", "Ciliated": "ciliated",
                       "Transitional Ciliated": "ciliated"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def composition_check(ad, arm_col: str) -> pd.DataFrame:
    """Called vs published per-sample composition, on the three coarse classes
    the published annotation can express (club and goblet are pooled, since the
    deposit does not split them)."""
    pub = G.published_composition(RAW / "GSE233145_cells_metadata.txt.gz")
    pub_coarse = (pub.T.groupby(PUBLISHED_TO_COARSE).sum().T)
    called = pd.crosstab(ad.obs["sample_id"], ad.obs[arm_col])
    called = called.div(called.sum(1), axis=0)
    called["secretory"] = called.reindex(columns=["club", "goblet"]).fillna(0).sum(1)
    rows = []
    for s in called.index:
        if s not in pub_coarse.index:
            continue
        for c in ["basal", "secretory", "ciliated"]:
            rows.append({"sample_id": s, "class": c,
                         "published": round(float(pub_coarse.loc[s, c]), 4),
                         "called": round(float(called.loc[s].get(c, 0.0)), 4)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--health-state", default="healthy", choices=["healthy", "COPD", "all"])
    args = ap.parse_args()
    state = None if args.health_state == "all" else args.health_state
    tag = args.health_state

    h5 = OUT / ("ali.h5ad" if tag == "healthy" else f"ali_{tag}.h5ad")
    if h5.exists():
        ad = sc.read_h5ad(h5)
        log(f"cached object: {ad.n_obs} cells")
    else:
        st = G.sample_table(ROOT / "data/external/gse233145_samples.tsv")
        ad = G.load(RAW / "raw", st, health_state=state, log=log)
        sc.pp.highly_variable_genes(ad, n_top_genes=N_HVG, batch_key="donor")
        ad.write(h5)
        log(f"object: {ad.n_obs} cells x {ad.n_vars} genes, "
            f"donors {sorted(ad.obs.donor.unique())}, days {sorted(ad.obs.day.unique())}")

    hvg = ad.var_names[ad.var.highly_variable].tolist()
    marker_genes = sorted({g for v in MARKERS.values() for g in
                           (v if isinstance(v, list) else sum(v.values(), []))})
    genes = sorted(set(hvg) | set(marker_genes))

    rfile = OUT / "hlca_reference.npz"
    if rfile.exists():
        z = np.load(rfile, allow_pickle=True)
        X_ref, y_ref, shared = z["X"], z["y"], list(z["genes"])
    else:
        X_ref, y_ref, shared = A.hlca_reference(genes, HLCA, log)
        np.savez_compressed(rfile, X=X_ref, y=y_ref, genes=np.array(shared))
    log(f"reference: {X_ref.shape[0]} atlas cells x {X_ref.shape[1]} genes")

    finest, margin, best = A.transfer(ad, X_ref, y_ref, shared, log)
    ad.obs["atlas_finest"] = finest
    ad.obs["transfer_margin"] = margin
    ad.obs["transfer_conf"] = best
    l3map = A.finest_to_l3(OUT / "finest_to_l3.csv", HLCA)
    ad.obs["atlas_l3"] = finest.map(l3map).fillna("unknown").values
    ad.obs["class_atlas"] = A.coarse(finest, l3map).values

    # Drop-seq at this depth never captured BPIFB2/MUC19, so the SMG-mucous panel
    # is unscoreable -- and an ALI culture grown from expanded basal cells has no
    # submucosal gland cells to veto in the first place. Declared explicitly so
    # the resolver records the drop instead of silently scoring a depleted panel.
    params = {**PARAMS, "secretory": {**PARAMS["secretory"], "optional_programs": ["smg"]}}
    sec = ad[ad.obs.atlas_l3.values == "Secretory"].copy()
    assign, fitted = R.resolve_secretory(sec, MARKERS, params, dataset_key="donor")
    ad.obs["class_resolver"] = ad.obs["class_atlas"].astype(object)
    ad.obs.loc[assign.index, "class_resolver"] = assign["cell_class"].values

    cols = ["sample", "sample_id", "donor", "day", "health_state", "atlas_finest",
            "atlas_l3", "transfer_margin", "transfer_conf", "class_atlas", "class_resolver"]
    ad.obs[cols].to_csv(OUT / ("cell_classes.csv" if tag == "healthy" else f"cell_classes_{tag}.csv"))
    ad.write(h5)

    summary = {"n_cells": int(ad.n_obs), "health_state": tag,
               "resolver_thresholds": fitted,
               "atlas_finest_counts": ad.obs.atlas_finest.value_counts().to_dict()}
    for arm in ["class_atlas", "class_resolver"]:
        tab = pd.crosstab(ad.obs["sample"], ad.obs[arm])
        tab.to_csv(OUT / f"composition_{arm}.csv")
        log(f"\n{arm}\n{tab.to_string()}")
        mp = A.marker_positivity(ad, ad.obs[arm])
        mp.to_csv(OUT / f"marker_positivity_{arm}.csv", index=False)
        log(f"\nmarker positivity ({arm})\n{mp.to_string(index=False)}")
        cc = composition_check(ad, arm)
        cc.to_csv(OUT / f"published_vs_called_composition_{arm}.csv", index=False)
        r = cc.groupby("class").apply(lambda d: d.published.corr(d.called), include_groups=False)
        mae = cc.groupby("class").apply(lambda d: (d.published - d.called).abs().mean(),
                                        include_groups=False)
        summary[f"composition_vs_published_{arm}"] = {
            "pearson_r": r.round(3).to_dict(), "mean_abs_diff": mae.round(3).to_dict()}
        log(f"\nvs published composition ({arm}) r={r.round(2).to_dict()} "
            f"MAE={mae.round(3).to_dict()}")
    json.dump(summary, open(OUT / f"prep_summary_{tag}.json", "w"), indent=1, default=float)
    log("done")


if __name__ == "__main__":
    main()
