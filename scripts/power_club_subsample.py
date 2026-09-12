"""Sample-size analysis: how reproducible is a VIPER TF ranking as a function
of the number of cells available in the case arm?

Design (split-half, so no ground truth is needed):
  * The control arm (basal) is held FULL and identical in both halves, so the
    only varying quantity is the number of case-lineage cells.
  * For each target N, case cells are drawn WITHOUT replacement into two
    DISJOINT halves of N cells each, allocated across (donor, day) groups in
    proportion to the cells each group has, so both halves keep the real
    donor/day composition of the experiment.
  * Each half goes through the production path unchanged: pseudobulk by
    (donor, day, class) with min_cells, limma-voom case-vs-basal blocked by
    donor, rank-to-normal-quantile GES, aREA against the same regulon.
  * Reproducibility = Spearman rho and top-20 overlap BETWEEN the two halves.
    This is the quantity an experiment must control: two independent samples
    of the same size from the same biology should agree on their top TFs.
  * The ciliated lineage is run as a positive control: its driver (FOXJ1) is
    known, so its curve shows what an adequately-powered contrast looks like
    in the same dataset with the same machinery.
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

os.environ.setdefault("NUMBA_CACHE_DIR", str(Path.cwd() / ".numba_cache"))

from basal_to_club.signature.pseudobulk import (make_pseudobulk, limma_voom_signature, zscore_signature)
from basal_to_club.viper.run_viper import run_pyviper

MIN_CELLS = 25          # pseudobulk group floor, as in the production run
MIN_SAMPLES = 2         # a contrast arm needs >= 2 pseudobulk samples
N_GRID = [50, 100, 200, 400, 800, 1500]
N_REPS = 3
CONTROLS = {"club": ["ELF3", "SOX2", "KLF5", "CEBPB", "FOXA2"],
            "ciliated": ["FOXJ1", "RFX2", "RFX3", "MYB", "TP73"]}


def load_regulon(tfgene: Path) -> pd.DataFrame:
    reg = pd.read_csv(tfgene, sep="\t")
    reg.columns = [c.lower() for c in reg.columns]
    reg = reg.rename(columns={"gene": "target", "score": "likelihood"})
    reg["mo"] = 1.0
    return reg[["tf", "target", "mo", "likelihood"]]


def allocate(sizes: pd.Series, n: int, rng) -> dict:
    """Split n cells across groups in proportion to size, capped by size."""
    frac = sizes / sizes.sum()
    take = np.floor(frac * n).astype(int)
    short = n - take.sum()
    if short > 0:                     # hand the remainder to the largest groups
        order = (frac * n - take).sort_values(ascending=False).index
        for g in order[:short]:
            take[g] += 1
    return {g: int(min(v, sizes[g])) for g, v in take.items()}


def rank_one(adata, case_cells, base_cells, lineage, regulon, label_col):
    """Full production path on one cell subset -> Series of rank percentiles."""
    sub = adata[list(case_cells) + list(base_cells)]
    counts, meta = make_pseudobulk(sub, ["donor", "day", label_col],
                                   min_cells=MIN_CELLS, layer="counts")
    meta = meta.rename(columns={label_col: "cell_class"})
    n_case = (meta.cell_class == lineage).sum()
    if n_case < MIN_SAMPLES:
        return None, n_case
    res = limma_voom_signature(counts, meta, "cell_class", lineage, "basal", block="donor")
    ges = zscore_signature(res["t_stat"]).to_frame(name="sig")
    act = run_pyviper(ges, regulon, n_perm=1000, pleiotropy=False)["sig"].dropna()
    pct = act.rank(ascending=False, method="min") / len(act)
    return pct, n_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drem-set", default="drem_ali")
    ap.add_argument("--label-col", default="class_atlas")
    ap.add_argument("--lineages", default="club,ciliated")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    import anndata as ad
    root = Path("results") / a.drem_set
    adata = ad.read_h5ad(root / "ali.h5ad")
    regulon = load_regulon(root / "tfgene_coexpression.txt")
    out = Path(a.out or f"results/power/{a.drem_set}")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    obs = adata.obs
    base_cells = obs.index[obs[a.label_col] == "basal"]
    rows, ref_rows = [], []

    for lineage in a.lineages.split(","):
        cells = obs.index[obs[a.label_col] == lineage]
        grp = obs.loc[cells].groupby(["donor", "day"], observed=True).size()
        grp = grp[grp > 0]
        total = int(grp.sum())
        ref, n_ref = rank_one(adata, cells, base_cells, lineage, regulon, a.label_col)
        if ref is not None:
            ref_rows += [{"lineage": lineage, "tf": t, "pct": p, "n_cells": total,
                          "n_samples": int(n_ref)} for t, p in ref.items()]
        print(f"[{lineage}] total={total} groups={len(grp)} ref_samples={n_ref}", flush=True)

        for n in [x for x in N_GRID if 2 * x <= total] + ([total // 2] if total // 2 not in N_GRID else []):
            for rep in range(N_REPS):
                alloc = allocate(grp, 2 * n, rng)
                half_a, half_b = [], []
                for (d, day), k in alloc.items():
                    pool = obs.index[(obs[a.label_col] == lineage) & (obs.donor == d) & (obs.day == day)]
                    pick = rng.choice(pool, size=k, replace=False)
                    half_a += list(pick[: k // 2]); half_b += list(pick[k // 2: 2 * (k // 2)])
                pa, na = rank_one(adata, half_a, base_cells, lineage, regulon, a.label_col)
                pb, nb = rank_one(adata, half_b, base_cells, lineage, regulon, a.label_col)
                rec = {"lineage": lineage, "n_cells": n, "rep": rep,
                       "n_samples_a": int(na), "n_samples_b": int(nb),
                       "ok": pa is not None and pb is not None}
                if rec["ok"]:
                    common = pa.index.intersection(pb.index)
                    rec["n_tf"] = len(common)
                    rec["spearman"] = float(pa[common].corr(pb[common], method="spearman"))
                    rec["top20_overlap"] = len(set(pa.sort_values().index[:20]) & set(pb.sort_values().index[:20]))
                    for tf in CONTROLS[lineage]:
                        if tf in pa.index and tf in pb.index:
                            rec[f"pct_{tf}_a"], rec[f"pct_{tf}_b"] = float(pa[tf]), float(pb[tf])
                    if ref is not None:
                        c2 = pa.index.intersection(ref.index)
                        rec["spearman_vs_full_a"] = float(pa[c2].corr(ref[c2], method="spearman"))
                        rec["top20_vs_full_a"] = len(set(pa.sort_values().index[:20]) & set(ref.sort_values().index[:20]))
                rows.append(rec)
                print(f"  n={n} rep={rep} ok={rec['ok']} rho={rec.get('spearman')} top20={rec.get('top20_overlap')}", flush=True)
                pd.DataFrame(rows).to_csv(out / "power_curve.csv", index=False)
    pd.DataFrame(ref_rows).to_csv(out / "reference_rankings.csv", index=False)
    json.dump({"drem_set": a.drem_set, "label_col": a.label_col, "min_cells": MIN_CELLS,
               "n_reps": N_REPS, "n_grid": N_GRID, "regulon_edges": int(len(regulon)),
               "regulon_tfs": int(regulon.tf.nunique())},
              open(out / "power_config.json", "w"), indent=1)


if __name__ == "__main__":
    main()
