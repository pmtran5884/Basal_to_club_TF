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


def hlca_reference(genes: list[str]):
    """Sampled HLCA airway-epithelial cells (log-normalised) for classifier training.

    Centroid-argmax transfer was tried first and failed an elementary check: it
    called 1,515 goblet cells in a sample with 0.5% MUC5AC+ cells. A regularised
    multinomial classifier on the same reference is used instead.
    """
    ref = sc.read_h5ad(HLCA, backed="r")
    obs = ref.obs
    name_col = "feature_name" if "feature_name" in ref.var.columns else None
    var_names = (ref.var[name_col].astype(str).values if name_col
                 else np.asarray(ref.var_names, dtype=str))
    gidx = pd.Series(np.arange(len(var_names)), index=var_names)
    gidx = gidx[~gidx.index.duplicated()]
    shared = [g for g in genes if g in gidx.index]
    cols = gidx.loc[shared].to_numpy()
    log(f"HLCA: {ref.n_obs} cells, {len(shared)}/{len(genes)} genes shared")

    keep = obs["ann_level_3"].astype(str).isin(AIRWAY_L3)
    rng = np.random.default_rng(0)
    blocks, labels = [], []
    for lab, idx in obs.index[keep].to_series().groupby(
            obs.loc[keep, "ann_finest_level"].astype(str), observed=True).groups.items():
        pos = np.flatnonzero(obs.index.isin(idx))
        if len(pos) < 100:
            continue
        if len(pos) > N_REF_PER_CLASS:
            pos = np.sort(rng.choice(pos, N_REF_PER_CLASS, replace=False))
        X = ref.X[pos, :]
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        blocks.append(X[:, cols].astype(np.float32))
        labels += [lab] * len(pos)
        log(f"  reference {lab}: {len(pos)} cells")
    ref.file.close()
    return np.vstack(blocks), np.array(labels), shared


def transfer(ad, X_ref, y_ref, genes):
    """L2 multinomial logistic regression trained on the atlas, applied to ALI."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split

    use = [i for i, g in enumerate(genes) if g in set(ad.var_names)]
    gnames = [genes[i] for i in use]
    X_ref = X_ref[:, use]
    mu, sd = X_ref.mean(0), X_ref.std(0) + 1e-9
    Z = (X_ref - mu) / sd
    Xtr, Xte, ytr, yte = train_test_split(Z, y_ref, test_size=0.2, stratify=y_ref,
                                          random_state=0)
    clf = LogisticRegression(C=0.1, max_iter=2000)
    clf.fit(Xtr, ytr)
    log(f"classifier: {X_ref.shape[0]} atlas cells x {len(gnames)} genes, "
        f"held-out accuracy {clf.score(Xte, yte):.3f}")
    clf.fit(Z, y_ref)
    Xq = ad[:, gnames].X
    Xq = np.asarray(Xq.todense()) if hasattr(Xq, "todense") else np.asarray(Xq)
    P = clf.predict_proba((Xq - mu) / sd)
    lab = pd.Series(clf.classes_[P.argmax(1)], index=ad.obs_names)
    conf = pd.Series(P.max(1), index=ad.obs_names)
    margin = pd.Series(np.sort(P, 1)[:, -1] - np.sort(P, 1)[:, -2], index=ad.obs_names)
    return lab, margin, conf


def marker_positivity(ad, cls: pd.Series) -> pd.DataFrame:
    """Fraction of cells in each called class expressing the canonical marker of
    that class and of its neighbours -- the check the centroid transfer failed."""
    probe = {"basal": "KRT5", "club": "SCGB1A1", "goblet": "MUC5AC",
             "goblet_alt": "MUC5B", "ciliated": "FOXJ1"}
    rows = []
    for c in sorted(cls.dropna().unique()):
        idx = cls.index[cls == c]
        sub = ad[idx]
        r = {"class": c, "n_cells": len(idx)}
        for k, g in probe.items():
            if g in ad.var_names:
                x = sub[:, g].X
                x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
                r[f"pct_{g}+"] = round(float((x > 0).mean()) * 100, 1)
        rows.append(r)
    return pd.DataFrame(rows)


def finest_to_l3() -> pd.Series:
    """ann_finest_level -> ann_level_3 from the atlas itself, so the coarse
    classes match results/atlas_label_arms exactly."""
    f = OUT / "finest_to_l3.csv"
    if f.exists():
        return pd.read_csv(f, index_col=0).iloc[:, 0]
    ref = sc.read_h5ad(HLCA, backed="r")
    m = (ref.obs[["ann_finest_level", "ann_level_3"]].astype(str)
            .drop_duplicates().set_index("ann_finest_level")["ann_level_3"])
    ref.file.close()
    m.to_csv(f)
    return m


def coarse(finest: pd.Series, l3map: pd.Series) -> pd.Series:
    """Same rule as the VIPER arms run: basal and ciliated from ann_level_3,
    club/goblet from the finest label. Cultured basal cells transfer mostly to
    Hillock-like, which is ann_level_3 'Basal'."""
    l3 = finest.map(l3map).fillna("unknown")
    cls = pd.Series("other", index=finest.index, dtype=object)
    cls[l3 == "Basal"] = "basal"
    cls[l3 == "Multiciliated lineage"] = "ciliated"
    cls[finest.isin(CLUB_LABELS)] = "club"
    cls[finest.isin(GOBLET_LABELS)] = "goblet"
    return cls


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
        X_ref, y_ref, shared = hlca_reference(genes)
        np.savez_compressed(rfile, X=X_ref, y=y_ref, genes=np.array(shared))
    log(f"reference: {X_ref.shape[0]} cells x {X_ref.shape[1]} genes, "
        f"{len(set(y_ref))} atlas labels")

    finest, margin, best = transfer(ad, X_ref, y_ref, shared)
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
