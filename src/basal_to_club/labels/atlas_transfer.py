"""Atlas label transfer and the coarse class mapping, shared by the dataset arms.

Lifted verbatim from scripts/prep_drem_ali.py so the nasal (GSE121600) and
bronchial (GSE233145) time courses are labelled by identical code. Centroid
-argmax transfer was tried first and failed an elementary check -- it called
1,515 goblet cells in a sample with 0.5% MUC5AC+ cells -- so a regularised
multinomial classifier trained on sampled HLCA core cells is used instead.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

AIRWAY_L3 = ("Basal", "Secretory", "Multiciliated lineage")
CLUB_LABELS = {"Club (non-nasal)", "Club (nasal)"}
GOBLET_LABELS = {"Goblet (nasal)", "Goblet (bronchial)", "Goblet (subsegmental)"}
N_REF_PER_CLASS = 2500


def hlca_reference(genes: list[str], hlca_path: Path, log=print):
    """Sampled HLCA airway-epithelial cells (log-normalised) for classifier training.

    Centroid-argmax transfer was tried first and failed an elementary check: it
    called 1,515 goblet cells in a sample with 0.5% MUC5AC+ cells. A regularised
    multinomial classifier on the same reference is used instead.
    """
    ref = sc.read_h5ad(hlca_path, backed="r")
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


def transfer(ad, X_ref, y_ref, genes, log=print):
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


def finest_to_l3(cache: Path, hlca_path: Path) -> pd.Series:
    """ann_finest_level -> ann_level_3 from the atlas itself, so the coarse
    classes match results/atlas_label_arms exactly."""
    f = Path(cache)
    if f.exists():
        return pd.read_csv(f, index_col=0).iloc[:, 0]
    ref = sc.read_h5ad(hlca_path, backed="r")
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
