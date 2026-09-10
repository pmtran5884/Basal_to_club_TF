"""Do the club/goblet class definitions change which TFs come out?

WHY THIS RUN EXISTS. `run_hlca_prior_only.py` defines the secretory classes with
the marker-score resolver. QC (docs/qc_2026-09/QC_REPORT.md) showed the resolver's
classes overlap the atlas's own finest-level labels only weakly -- goblet Jaccard
0.02-0.20 per study -- so the TF rankings it produced may be an artefact of that
particular class definition rather than of club vs goblet biology.

WHAT IT DOES. Two arms on the SAME cells, the SAME donors, the SAME studies and
the SAME regulons, differing only in how club and goblet are defined:

    arm "atlas_label"  club/goblet from HLCA ann_finest_level
    arm "resolver"     club/goblet from the marker-score resolver

Studies are restricted to those that delineate both classes with enough cells to
support a comparison (>= MIN_CELLS_PER_CLASS atlas-labelled cells in club AND in
goblet), so that a difference between arms cannot be attributed to study mix.

Regulons are NOT rebuilt: the ChEA3 networks cached by the prior-only run are
reused verbatim, which holds the network fixed and leaves the signature as the
only moving part. Consequence: the gate still cannot pass (no de novo family),
so every ranking here remains provisional, exactly as in the prior-only run.

Basal and ciliated are taken from the atlas label in BOTH arms (`ann_level_3`);
only the secretory compartment is in question.
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
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from basal_to_club.signature import pseudobulk as PB          # noqa: E402
from basal_to_club.viper.run_metaviper import integrate        # noqa: E402
from basal_to_club.viper.run_viper import fdr as nes_fdr       # noqa: E402
from basal_to_club.viper.run_viper import run_pyviper          # noqa: E402

PRIOR_RUN = ROOT / "results" / "hlca_prior_only"
RES = ROOT / "results" / "atlas_label_arms"
H5 = ROOT / "data/raw/hlca_core/hlca_core.h5ad"

PARAMS = yaml.safe_load(open(ROOT / "config/params.yaml"))
AIRWAY_L3 = ("Basal", "Secretory", "Multiciliated lineage", "Submucosal Secretory", "Rare")

CLUB_LABELS = {"Club (non-nasal)", "Club (nasal)"}
GOBLET_LABELS = {"Goblet (nasal)", "Goblet (bronchial)", "Goblet (subsegmental)"}
MIN_CELLS_PER_CLASS = 100

# ChEA3 library families. The gate's "both families" criterion means binding
# evidence AND expression evidence, not two libraries of the same kind.
BINDING = {"ENCODE_ChIP-seq", "ReMap_ChIP-seq", "Literature_ChIP-seq"}
CONTRASTS = ("club", "goblet", "ciliated")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_subset():
    """Same subsetting rule as the prior-only run, so the two are comparable."""
    cache = PRIOR_RUN / "airway_subset.h5ad"
    if cache.exists():
        log(f"loading cached subset {cache}")
        return sc.read_h5ad(cache)
    log("reading HLCA backed")
    ad = sc.read_h5ad(H5, backed="r")
    obs = ad.obs
    m = (obs["ann_level_3"].isin(AIRWAY_L3).to_numpy()
         & (obs["disease"] == "normal").to_numpy()
         & (obs["tissue"] != "nose").to_numpy())
    log(f"subsetting {int(m.sum())} cells of {ad.n_obs}")
    sub = ad[m].to_memory()
    cache.parent.mkdir(parents=True, exist_ok=True)
    sub.write_h5ad(cache)
    return sub


def prepare(ad):
    """CELLxGENE convention: .raw holds raw counts, X is normalized."""
    counts = ad.raw.to_adata() if ad.raw is not None else ad.copy()
    counts.var_names = ad.raw.var["feature_name"].astype(str).to_numpy() if ad.raw is not None \
        else ad.var["feature_name"].astype(str).to_numpy()
    counts.var_names_make_unique()
    counts.obs = ad.obs.copy()
    return counts


def trustworthy_studies(obs: pd.DataFrame) -> list[str]:
    fin = obs["ann_finest_level"].astype(str)
    tab = pd.DataFrame({
        "study": obs["study"].astype(str),
        "is_club": fin.isin(CLUB_LABELS),
        "is_goblet": fin.isin(GOBLET_LABELS),
    }).groupby("study")[["is_club", "is_goblet"]].sum()
    keep = tab[(tab.is_club >= MIN_CELLS_PER_CLASS) & (tab.is_goblet >= MIN_CELLS_PER_CLASS)]
    tab.to_csv(RES / "study_class_counts.csv")
    return sorted(keep.index)


def classes_atlas_label(obs: pd.DataFrame) -> pd.Series:
    l3 = obs["ann_level_3"].astype(str)
    fin = obs["ann_finest_level"].astype(str)
    cls = pd.Series("other", index=obs.index, dtype=object)
    cls[l3 == "Basal"] = "basal"
    cls[l3 == "Multiciliated lineage"] = "ciliated"
    cls[fin.isin(CLUB_LABELS)] = "club"
    cls[fin.isin(GOBLET_LABELS)] = "goblet"
    return cls


def classes_resolver(obs: pd.DataFrame) -> pd.Series:
    """Resolver assignments as cached by the prior-only run (not refitted)."""
    assign = pd.read_csv(PRIOR_RUN / "resolver_assignments.csv", index_col=0)
    l3 = obs["ann_level_3"].astype(str)
    cls = pd.Series("other", index=obs.index, dtype=object)
    cls[l3 == "Basal"] = "basal"
    cls[l3 == "Multiciliated lineage"] = "ciliated"
    shared = obs.index.intersection(assign.index)
    cls.loc[shared] = assign.loc[shared, "cell_class"].astype(str)
    return cls


def signatures_for(counts, cls: pd.Series, arm_dir: Path) -> pd.DataFrame:
    keep = cls.isin(["basal", *CONTRASTS])
    ad = counts[keep.to_numpy()].copy()
    ad.obs["cell_class"] = cls[keep].to_numpy()
    ad.obs["dataset_id"] = ad.obs["study"].astype(str)
    ad.layers["counts"] = ad.X
    pbc, meta = PB.make_pseudobulk(
        ad, ["dataset_id", "donor_id", "cell_class"],
        min_cells=PARAMS["signature"]["pseudobulk"]["min_cells_per_pseudobulk"])
    meta.to_csv(arm_dir / "pseudobulk_meta.csv")
    log(f"  pseudobulk profiles={pbc.shape[0]}\n"
        + meta["cell_class"].value_counts().to_string())
    sigs = {}
    for target in CONTRASTS:
        sub = meta["cell_class"].isin(["basal", target])
        res = PB.limma_voom_signature(pbc.loc[sub.values], meta.loc[sub.values],
                                      group_col="cell_class", case=target, control="basal")
        res.to_csv(arm_dir / f"limma_basal_to_{target}.csv")
        sigs[f"basal_to_{target}"] = PB.zscore_signature(res["t_stat"])
        log(f"  signature basal_to_{target}: {len(res)} genes")
    sig = pd.DataFrame(sigs)
    sig.to_csv(arm_dir / "signatures_zscore.csv")
    return sig


def viper_for(sig: pd.DataFrame, arm_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score every cached regulon, then integrate across networks (metaVIPER)."""
    per_lib = []
    for path in sorted(PRIOR_RUN.glob("regulon_*.csv")):
        name = path.stem.replace("regulon_", "")
        reg = pd.read_csv(path)
        act = run_pyviper(sig, reg[["tf", "target", "mo", "likelihood"]],
                          n_perm=0, pleiotropy=False)
        act.to_csv(arm_dir / f"viper_nes_{name}.csv")
        m = (act.reset_index().rename(columns={"index": "tf"})
                .melt(id_vars="tf", var_name="contrast", value_name="nes")
                .dropna(subset=["nes"]))
        m["network"], m["method"] = name, "viper"
        # per-network FDR, as the prior-only run does; integrate() takes its min
        m["fdr"] = m.groupby("contrast")["nes"].transform(nes_fdr)
        per_lib.append(m)
        log(f"  viper {name}: {act.shape[0]} TFs")
    long = pd.concat(per_lib, ignore_index=True)
    long.to_csv(arm_dir / "tf_activity_long.csv", index=False)
    comb = integrate(long, ["tf", "contrast", "method"])
    comb.to_csv(arm_dir / "metaviper_combined.csv", index=False)
    return long, comb


def rank_table(long: pd.DataFrame, comb: pd.DataFrame) -> pd.DataFrame:
    """Rank per contrast, and flag two-family support as the gate defines it."""
    fam = long.assign(family=np.where(long.network.isin(BINDING), "binding", "expression"))
    pos = (fam.assign(pos=fam.nes > 0)
              .groupby(["tf", "contrast", "family"])["pos"].any().unstack("family"))
    for c in ("binding", "expression"):
        if c not in pos:
            pos[c] = False
    nlib = long.dropna(subset=["nes"]).groupby(["tf", "contrast"])["network"].nunique()
    out = comb.set_index(["tf", "contrast"]).join(
        pos.fillna(False).rename(columns=lambda c: f"{c}_positive")).join(nlib.rename("n_libraries"))
    out["both_families_positive"] = out.binding_positive & out.expression_positive
    out = out.reset_index()
    out["rank"] = out.groupby("contrast")["nes"].rank(ascending=False).astype(int)
    return out.sort_values(["contrast", "rank"])


def main():
    RES.mkdir(parents=True, exist_ok=True)
    counts = prepare(load_subset())
    studies = trustworthy_studies(counts.obs)
    log(f"trustworthy studies (>= {MIN_CELLS_PER_CLASS} atlas cells in both classes): {studies}")
    counts = counts[counts.obs["study"].astype(str).isin(studies)].copy()
    log(f"cells retained {counts.n_obs} genes {counts.n_vars}")

    arms = {"atlas_label": classes_atlas_label(counts.obs),
            "resolver": classes_resolver(counts.obs)}
    ranks = {}
    for arm, cls in arms.items():
        arm_dir = RES / arm
        arm_dir.mkdir(parents=True, exist_ok=True)
        log(f"=== arm {arm} ===\n" + cls.value_counts().to_string())
        cls.to_frame("cell_class").to_csv(arm_dir / "cell_classes.csv")
        sig = signatures_for(counts, cls, arm_dir)
        long, comb = viper_for(sig, arm_dir)
        rk = rank_table(long, comb)
        rk.to_csv(arm_dir / "tf_ranking.csv", index=False)
        ranks[arm] = rk

    # ---- did the definition change the answer? -------------------------------
    rows = []
    for contrast in [f"basal_to_{c}" for c in CONTRASTS]:
        a = ranks["atlas_label"].query("contrast == @contrast").set_index("tf")["nes"]
        b = ranks["resolver"].query("contrast == @contrast").set_index("tf")["nes"]
        shared = a.index.intersection(b.index)
        ta, tb = set(a.nlargest(20).index), set(b.nlargest(20).index)
        rows.append({"contrast": contrast, "n_tf_shared": len(shared),
                     "spearman_rho": spearmanr(a[shared], b[shared]).statistic,
                     "top20_overlap": len(ta & tb),
                     "new_in_atlas_label_top20": ",".join(sorted(ta - tb)),
                     "lost_from_resolver_top20": ",".join(sorted(tb - ta))})
    cmp = pd.DataFrame(rows)
    cmp.to_csv(RES / "arm_comparison.csv", index=False)
    log("arm comparison\n" + cmp.drop(columns=["new_in_atlas_label_top20",
                                               "lost_from_resolver_top20"]).to_string(index=False))

    # ---- club-specific vs the goblet contrast, per arm -----------------------
    uniq = {}
    for arm, rk in ranks.items():
        club = rk.query("contrast == 'basal_to_club'").nlargest(20, "nes")
        gob = set(rk.query("contrast == 'basal_to_goblet'").nlargest(20, "nes")["tf"])
        club["in_goblet_top20"] = club.tf.isin(gob)
        club.to_csv(RES / arm / "club_top20_vs_goblet.csv", index=False)
        uniq[arm] = club.loc[~club.in_goblet_top20, "tf"].tolist()
        log(f"{arm}: club top-20 not in goblet top-20 -> {uniq[arm]}")
    (RES / "club_unique_tfs.json").write_text(json.dumps(uniq, indent=1))
    log("done")


if __name__ == "__main__":
    main()
