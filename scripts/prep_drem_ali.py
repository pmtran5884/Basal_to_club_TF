"""Stage 1 for the DREM arm: build the ALI time-course object and both label arms.

Data: Ruiz Garcia et al. 2019 (GSE121600), nasal ALI cultures, the only planned
dataset with a true *differentiation* time course. Donor D246 days 2-22,
donor D275 days 7-47. The deposit carries raw droplet matrices only, so both
label arms have to be produced here.

Arms (design mirrors results/atlas_label_arms so DREM and VIPER are compared
under the same two class definitions):
  atlas_label : HLCA core labels transferred by centroid correlation; the
                atlas's own Club/Goblet split is kept.
  resolver    : basal and ciliated identical to the atlas arm; the secretory
                compartment is re-split by the club/goblet resolver.
"""
from __future__ import annotations

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
from basal_to_club.secretory import resolver as R  # noqa: E402
from basal_to_club.labels.atlas_transfer import (  # noqa: E402
    AIRWAY_L3, CLUB_LABELS, GOBLET_LABELS, coarse, hlca_reference, marker_positivity, transfer,
    finest_to_l3 as _finest_to_l3)


def finest_to_l3():
    return _finest_to_l3(OUT / "finest_to_l3.csv", HLCA)

RAW = ROOT / "data/raw/ali_gse121600"
HLCA = ROOT / "data/raw/hlca_core/hlca_core.h5ad"
OUT = ROOT / "results/drem_ali"
OUT.mkdir(parents=True, exist_ok=True)

MARKERS = yaml.safe_load(open(ROOT / "config/markers.yaml"))["markers"]
PARAMS = yaml.safe_load(open(ROOT / "config/params.yaml"))

AIRWAY_L3 = ("Basal", "Secretory", "Multiciliated lineage")
CLUB_LABELS = {"Club (non-nasal)", "Club (nasal)"}
GOBLET_LABELS = {"Goblet (nasal)", "Goblet (bronchial)", "Goblet (subsegmental)"}
N_REF_PER_CLASS = 2500
N_HVG = 3000


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_ali() -> "sc.AnnData":
    """Per-sample QC then concatenation. Sample names encode donor and ALI day."""
    ads = []
    for f in sorted(RAW.glob("*.h5")):
        stem = f.name.split("_")[1].replace(".h5", "")      # e.g. D246ali2
        if "ali" not in stem.lower():
            log(f"skip {stem}: not an ALI time point")
            continue
        donor, day = stem.lower().split("ali")
        a = sc.read_10x_h5(f)
        a.var_names_make_unique()
        sc.pp.filter_cells(a, min_counts=1000)
        sc.pp.filter_cells(a, min_genes=500)
        mt = [g for g in a.var_names if g.startswith("MT-")]
        pmt = np.asarray(a[:, mt].X.sum(1)).ravel() / np.asarray(a.X.sum(1)).ravel() * 100
        a = a[pmt < 25].copy()
        a.obs["donor"] = donor.upper()
        a.obs["day"] = int(day)
        a.obs["sample"] = stem
        a.obs_names = [f"{stem}_{b}" for b in a.obs_names]
        ads.append(a)
        log(f"{stem}: {a.n_obs} cells")
    ad = sc.concat(ads, join="inner")
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad


def main():
    h5 = OUT / "ali.h5ad"
    if h5.exists():
        ad = sc.read_h5ad(h5)
        log(f"cached ALI object: {ad.n_obs} cells")
    else:
        ad = load_ali()
        sc.pp.highly_variable_genes(ad, n_top_genes=N_HVG, batch_key="donor")
        ad.write(h5)
        log(f"ALI object: {ad.n_obs} cells x {ad.n_vars} genes")

    hvg = ad.var_names[ad.var.highly_variable].tolist()
    marker_genes = sorted({g for v in MARKERS.values() for g in
                           (v if isinstance(v, list) else sum(v.values(), []))})
    genes = sorted(set(hvg) | set(marker_genes))

    rfile = OUT / "hlca_reference.npz"
    if rfile.exists():
        z = np.load(rfile, allow_pickle=True)
        X_ref, y_ref, shared = z["X"], z["y"], list(z["genes"])
    else:
        X_ref, y_ref, shared = hlca_reference(genes, HLCA, log)
        np.savez_compressed(rfile, X=X_ref, y=y_ref, genes=np.array(shared))
    log(f"reference: {X_ref.shape[0]} cells x {X_ref.shape[1]} genes, "
        f"{len(set(y_ref))} atlas labels")

    finest, margin, best = transfer(ad, X_ref, y_ref, shared, log)
    ad.obs["atlas_finest"] = finest
    ad.obs["transfer_margin"] = margin
    ad.obs["transfer_conf"] = best
    l3map = finest_to_l3()
    ad.obs["atlas_l3"] = finest.map(l3map).fillna("unknown").values
    ad.obs["class_atlas"] = coarse(finest, l3map).values

    # the resolver re-splits the whole secretory compartment, as in the VIPER arms
    sec = ad[ad.obs.atlas_l3.values == "Secretory"].copy()
    assign, fitted = R.resolve_secretory(sec, MARKERS, PARAMS, dataset_key="donor")
    ad.obs["class_resolver"] = ad.obs["class_atlas"].astype(object)
    ad.obs.loc[assign.index, "class_resolver"] = assign["cell_class"].values

    cls = ad.obs[["sample", "donor", "day", "atlas_finest", "atlas_l3",
                  "transfer_margin", "transfer_conf", "class_atlas", "class_resolver"]]
    cls.to_csv(OUT / "cell_classes.csv")
    ad.write(h5)

    for arm in ["class_atlas", "class_resolver"]:
        tab = pd.crosstab(ad.obs["sample"], ad.obs[arm])
        tab.to_csv(OUT / f"composition_{arm}.csv")
        log(f"\n{arm}\n{tab.to_string()}")
        mp = marker_positivity(ad, ad.obs[arm])
        mp.to_csv(OUT / f"marker_positivity_{arm}.csv", index=False)
        log(f"\nmarker positivity ({arm})\n{mp.to_string(index=False)}")
    json.dump({"resolver_thresholds": fitted,
               "n_cells": int(ad.n_obs),
               "atlas_finest_counts": ad.obs.atlas_finest.value_counts().to_dict()},
              open(OUT / "prep_summary.json", "w"), indent=1, default=float)
    log("done")


if __name__ == "__main__":
    main()
