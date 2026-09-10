"""Scoped run: HLCA core, prior-network family only.

WHAT THIS IS. A diagnostic execution of the pipeline's scientific core on one
already-integrated atlas. It enters the DAG below integration and runs:

    resolver -> pseudobulk -> limma-voom signatures -> ChEA3 regulons
             -> metaVIPER -> benchmark scoring

WHAT THIS IS NOT. It is not the gated pipeline. The gate in config/benchmark.yaml
requires `require_both_network_families: true`, and this run has only the prior
family (no ARACNe jar is installed). It therefore CANNOT pass the gate, by design,
and no club prediction produced here is gate-backed. Reported for diagnosis only.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from basal_to_club.secretory import resolver as R          # noqa: E402
from basal_to_club.signature import pseudobulk as PB        # noqa: E402
from basal_to_club.network import prior as PRIOR            # noqa: E402
from basal_to_club.benchmark import score as SCORE          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "hlca_prior_only"
RES.mkdir(parents=True, exist_ok=True)
H5 = ROOT / "data/raw/hlca_core/hlca_core.h5ad"
CHEA = ROOT / "data/external/chea3"

MARKERS = yaml.safe_load(open(ROOT / "config/markers.yaml"))["markers"]
PARAMS = yaml.safe_load(open(ROOT / "config/params.yaml"))
BENCH = yaml.safe_load(open(ROOT / "config/benchmark.yaml"))

AIRWAY_L3 = ("Basal", "Secretory", "Multiciliated lineage", "Submucosal Secretory", "Rare")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_subset():
    cache = RES / "airway_subset.h5ad"
    if cache.exists():
        log(f"loading cached subset {cache.name}")
        return sc.read_h5ad(cache)
    log("reading HLCA backed")
    ad = sc.read_h5ad(H5, backed="r")
    obs = ad.obs
    m = (obs["ann_level_3"].isin(AIRWAY_L3).to_numpy()
         & (obs["disease"] == "normal").to_numpy()
         & (obs["tissue"] != "nose").to_numpy())
    log(f"subsetting {int(m.sum())} cells of {ad.n_obs}")
    sub = ad[m].to_memory()
    sub.write_h5ad(cache)
    return sub


def prepare(ad):
    """CELLxGENE convention: .raw holds raw counts, X is normalized."""
    counts = ad.raw.to_adata() if ad.raw is not None else ad.copy()
    counts.var_names = ad.raw.var["feature_name"].astype(str).to_numpy() if ad.raw is not None \
        else ad.var["feature_name"].astype(str).to_numpy()
    counts.var_names_make_unique()
    counts.obs = ad.obs.copy()
    ln = counts.copy()
    sc.pp.normalize_total(ln, target_sum=1e4)
    sc.pp.log1p(ln)
    return counts, ln


def assign_classes(ad_ln, counts):
    """Basal / ciliated from the atlas label; secretory compartment re-derived."""
    l3 = ad_ln.obs["ann_level_3"].astype(str)
    cls = pd.Series("other", index=ad_ln.obs_names, dtype=object)
    cls[l3 == "Basal"] = "basal"
    cls[l3 == "Multiciliated lineage"] = "ciliated"

    sec_mask = l3.isin(["Secretory", "Submucosal Secretory"]).to_numpy()
    sec = ad_ln[sec_mask].copy()
    sec.obs["dataset_id"] = sec.obs["study"].astype(str)
    log(f"resolver on {sec.n_obs} secretory-compartment cells "
        f"across {sec.obs['dataset_id'].nunique()} studies")
    assign, fitted = R.resolve_secretory(sec, MARKERS, PARAMS, dataset_key="dataset_id")
    cls[assign.index] = assign["cell_class"].astype(str)

    counts.obs["cell_class"] = cls.reindex(counts.obs_names).to_numpy()
    ad_ln.obs["cell_class"] = cls.to_numpy()
    (RES / "resolver_thresholds.json").write_text(json.dumps(fitted, indent=1, default=float))
    assign.join(sec.obs[["ann_finest_level", "dataset_id"]]).to_csv(RES / "resolver_assignments.csv")
    return assign, fitted


def main():
    ad = load_subset()
    counts, ln = prepare(ad)
    log(f"genes={counts.n_vars} cells={counts.n_obs}")
    assign, fitted = assign_classes(ln, counts)

    xt = pd.crosstab(assign["cell_class"], ln.obs.loc[assign.index, "ann_finest_level"])
    xt.to_csv(RES / "resolver_vs_atlas_label.csv")
    log("class counts:\n" + counts.obs["cell_class"].value_counts().to_string())

    keep = counts.obs["cell_class"].isin(["basal", "club", "goblet", "ciliated"])
    pb_ad = counts[keep.to_numpy()].copy()
    pb_ad.obs["dataset_id"] = pb_ad.obs["study"].astype(str)
    pb_ad.layers["counts"] = pb_ad.X
    pbc, meta = PB.make_pseudobulk(pb_ad, ["dataset_id", "donor_id", "cell_class"],
                                   min_cells=PARAMS["signature"]["pseudobulk"]["min_cells_per_pseudobulk"])
    log(f"pseudobulk profiles={pbc.shape[0]} genes={pbc.shape[1]}")
    meta.to_csv(RES / "pseudobulk_meta.csv")
    log("profiles per class:\n" + meta["cell_class"].value_counts().to_string())

    sigs = {}
    for target in ["club", "goblet", "ciliated"]:
        sub = meta["cell_class"].isin(["basal", target])
        res = PB.limma_voom_signature(pbc.loc[sub.values], meta.loc[sub.values],
                                      group_col="cell_class", case=target, control="basal")
        sigs[f"basal_to_{target}"] = PB.zscore_signature(res["t_stat"])
        res.to_csv(RES / f"limma_basal_to_{target}.csv")
        log(f"signature basal_to_{target}: {len(res)} genes tested")
    sig_df = pd.DataFrame(sigs)
    sig_df.to_csv(RES / "signatures_zscore.csv")

    # ---- prior regulons from ChEA3 -------------------------------------------
    cpm = np.log1p(pbc.div(pbc.sum(axis=1), axis=0) * 1e6)
    libs = {p.stem: PRIOR.read_gmt(p) for p in sorted(CHEA.glob("*.gmt"))}
    support = PRIOR.cross_library_support(libs)
    regulons = {}
    for name, gs in libs.items():
        reg = PRIOR.build_regulon(gs, cpm, library_support=support,
                                  min_targets=25, max_targets=500)
        regulons[name] = reg
        log(f"regulon {name}: {reg['tf'].nunique()} TFs, {len(reg)} edges")
        reg.to_csv(RES / f"regulon_{name}.csv", index=False)

    # ---- VIPER ---------------------------------------------------------------
    from basal_to_club.viper.run_metaviper import integrate
    from basal_to_club.viper.run_viper import fdr as nes_fdr
    from basal_to_club.viper.run_viper import run_pyviper

    frames = []
    for name, reg in regulons.items():
        act = run_pyviper(sig_df, reg[["tf", "target", "mo", "likelihood"]],
                          n_perm=0, pleiotropy=False)
        act.to_csv(RES / f"viper_nes_{name}.csv")
        m = act.reset_index().rename(columns={"index": "tf"}).melt(
            id_vars="tf", var_name="contrast", value_name="nes").dropna(subset=["nes"])
        m["network"], m["method"], m["family"] = name, "viper", "prior"
        m["fdr"] = m.groupby("contrast")["nes"].transform(nes_fdr)
        frames.append(m)
        log(f"viper {name}: {act.shape[0]} TFs x {act.notna().sum().to_dict()}")

    activity = pd.concat(frames, ignore_index=True)
    activity.to_csv(RES / "tf_activity_long.csv", index=False)
    # metaVIPER integrates in LONG form: a TF present in 3 of 6 libraries is
    # combined over those 3. Intersecting the libraries instead would discard
    # every TF absent from the smallest ChIP-seq library, SPDEF included.
    comb = integrate(activity, ["tf", "contrast", "method"])
    comb.to_csv(RES / "metaviper_combined.csv", index=False)
    by_fam = integrate(activity, ["tf", "contrast", "method", "network"])
    by_fam.to_csv(RES / "metaviper_by_network.csv", index=False)
    log(f"metaVIPER: {comb.tf.nunique()} TFs across {comb.contrast.nunique()} contrasts")

    log("done")


if __name__ == "__main__":
    main()
