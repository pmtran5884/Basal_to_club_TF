"""Stage 2: DREM over the ALI differentiation time course.

One DREM model per (label arm x donor x lineage x regulon family). A lineage
series is the pseudobulk trajectory of basal cells plus that lineage only, so
the temporal signal is the emergence of that program out of the basal pool;
the whole-culture series (lineage "all") is the label-free reference model.

TF ranking: for each model, a TF's score is its most significant enrichment
p-value ("Score Split") on any *upward* edge, i.e. among genes induced along a
branch. Lower is stronger. This is the DREM analogue of a positive VIPER NES.
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from basal_to_club.drem import io as DIO      # noqa: E402
from basal_to_club.drem import parse as DP    # noqa: E402
from basal_to_club.drem import run as DR      # noqa: E402
from basal_to_club.network import prior as PRIOR  # noqa: E402

# Which time course to model: "drem_ali" (nasal, GSE121600) or "drem_bronchial"
# (GSE233145). Both prep scripts write ali.h5ad + cell_classes.csv into their
# own results dir, so the grid code is identical for the two datasets.
OUT = ROOT / "results" / os.environ.get("BTC_DREM_SET", "drem_ali")
CHEA = ROOT / "data/external/chea3"
JAR = ROOT.parent / "external/STEM_DREM/drem.jar"
CLASSES = ROOT.parent / "external/STEM_DREM/patched"
N_WORKERS = 5

FAMILIES = {
    "binding": ["ENCODE_ChIP-seq", "ReMap_ChIP-seq", "Literature_ChIP-seq"],
    "coexpression": ["ARCHS4_Coexpression", "GTEx_Coexpression", "Enrichr_Queries"],
}
ARMS = {"atlas_label": "class_atlas", "resolver": "class_resolver"}
LINEAGES = ["club", "goblet", "ciliated", "all"]
MIN_CELLS_PER_DAY = 20
# DREM's path search grows steeply with the number of sampled time points, so a
# dense course can be thinned to a representative grid via BTC_DREM_DAYS.
DAY_SUBSET = ({int(x) for x in os.environ["BTC_DREM_DAYS"].split(",")}
              if os.environ.get("BTC_DREM_DAYS") else None)
MIN_DAYS = 3
MAX_TARGETS_PER_TF = 1000
MIN_TARGETS_PER_TF = 15
MIN_CPM = 10.0


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def build_series(ad, cls: pd.Series, donor: str, lineage: str):
    """genes x days log2 CPM pseudobulk, plus the per-day cell counts."""
    keep = ad.obs.donor.values == donor
    if lineage != "all":
        keep &= cls.isin(["basal", lineage]).values
    sub = ad[keep]
    rows, counts = {}, []
    for day in sorted(sub.obs.day.unique()):
        if DAY_SUBSET is not None and int(day) not in DAY_SUBSET:
            continue
        cell = sub.obs.day.values == day
        n_lin = int((cls.loc[sub.obs_names[cell]] == lineage).sum()) if lineage != "all" else int(cell.sum())
        if cell.sum() < MIN_CELLS_PER_DAY:
            continue
        c = np.asarray(sub[cell].layers["counts"].sum(0)).ravel()
        rows[f"day{day}"] = np.log2(c / c.sum() * 1e6 + 1.0)
        counts.append({"donor": donor, "lineage": lineage, "day": int(day),
                       "n_cells": int(cell.sum()), "n_lineage": n_lin})
    if len(rows) < MIN_DAYS:
        return None, pd.DataFrame(counts)
    ser = pd.DataFrame(rows, index=ad.var_names)
    expressed = (2 ** ser - 1).max(1) >= MIN_CPM
    ser = ser[expressed]
    return ser, pd.DataFrame(counts)


def build_tf_gene(family: str, universe: set[str]) -> pd.DataFrame:
    """Union of the family's ChEA3 libraries, restricted to the expressed
    universe, capped per TF by cross-library support (ties broken by symbol)."""
    libs = {n: PRIOR.read_gmt(CHEA / f"{n}.gmt") for n in FAMILIES[family]}
    pairs: dict[tuple[str, str], int] = {}
    for gs in libs.values():
        for key, targets in gs.items():
            tf = key.split("_")[0].strip().upper()
            for g in targets:
                if g in universe:
                    pairs[(tf, g)] = pairs.get((tf, g), 0) + 1
    df = pd.DataFrame([{"tf": t, "target": g, "support": s} for (t, g), s in pairs.items()])
    df = (df.sort_values(["tf", "support", "target"], ascending=[True, False, True])
            .groupby("tf", group_keys=False).head(MAX_TARGETS_PER_TF))
    n = df.groupby("tf").size()
    df = df[df.tf.isin(n[n >= MIN_TARGETS_PER_TF].index)]
    # a TF that is not itself expressed in the culture cannot be driving it,
    # and every extra TF is a covariate in DREM's per-split regression
    df = df[df.tf.isin(universe)]
    return df


def _finished(d: Path) -> bool:
    """A model directory is complete when DREM printed its runtime and wrote
    at least one regulator table."""
    so = d / "drem_stdout.txt"
    if not so.exists() or "Time:" not in so.read_text():
        return False
    return any((d / "tfscores").glob("path_*.txt")) or any((d / "tfscores").glob("split_*.txt"))


def _one_run(task: dict) -> dict:
    """Execute one DREM model and rank its TFs.

    Runs in a worker thread: the cost is the DREM subprocess, and the sandbox
    does not permit the semaphores ProcessPoolExecutor needs.
    """
    d = Path(task["dir"])
    t0 = time.time()
    if _finished(d):
        proc = None
    else:
        proc = DR.run_drem(JAR, task["settings"], d / "model.txt", d / "geneassign.txt",
                           d / "tfscores", workdir=d, heap="4g", timeout=10800,
                           classes_dir=CLASSES)
        (d / "drem_stdout.txt").write_text((proc.stdout or "") + (proc.stderr or ""))
    ok = _finished(d)
    edges = DP.parse_edges(d / "tfscores")
    if not edges.empty:
        rk, source = DP.rank_tfs(edges), "path_tables"
    else:
        splits = DP.parse_splits(d / "tfscores")
        rk, source = DP.rank_tfs_from_splits(splits), "split_tables"
    if rk is not None and not rk.empty:
        rk.to_csv(d / "tf_ranking.csv", index=False)
    return {**{k: task[k] for k in ("arm", "donor", "lineage", "family",
                                    "n_genes", "n_days", "days")},
            "ok": bool(ok), "source": source, "refit": proc is not None,
            "n_tfs": int(0 if rk is None or rk.empty else rk.tf.nunique()),
            "seconds": round(time.time() - t0, 1), "_ranking": rk}


def main():
    ad = sc.read_h5ad(OUT / "ali.h5ad")
    cls_tbl = pd.read_csv(OUT / "cell_classes.csv", index_col=0)
    log(f"{ad.n_obs} cells; donors {sorted(ad.obs.donor.unique())}")
    log(f"day grid: {sorted(DAY_SUBSET) if DAY_SUBSET else 'all sampled days'}")

    # ---- series first, so the regulon universe is the expressed gene set ----
    series, counts = {}, []
    for arm, col in ARMS.items():
        cls = cls_tbl[col].reindex(ad.obs_names)
        for donor in sorted(ad.obs.donor.unique()):
            for lineage in LINEAGES:
                s, c = build_series(ad, cls, donor, lineage)
                counts.append(c)
                if s is None:
                    log(f"skip {arm}/{donor}/{lineage}: <{MIN_DAYS} usable days")
                    continue
                series[(arm, donor, lineage)] = s
    pd.concat(counts).drop_duplicates().to_csv(OUT / "series_cell_counts.csv", index=False)
    universe = set().union(*[set(s.index) for s in series.values()])
    log(f"{len(series)} series; {len(universe)} expressed genes")

    tfg = {}
    for fam in FAMILIES:
        p = OUT / f"tfgene_{fam}.txt"
        df = build_tf_gene(fam, universe)
        DIO.write_tf_gene(df, p)
        tfg[fam] = p
        log(f"{fam}: {df.tf.nunique()} TFs, {len(df)} edges -> {p.name}")

    # ---- one DREM model per series x family, run in parallel ----
    tasks = []
    for (arm, donor, lineage), s in sorted(series.items()):
        for fam, tf_file in tfg.items():
            d = OUT / "models" / arm / donor / lineage / fam
            d.mkdir(parents=True, exist_ok=True)
            if not _finished(d):
                expr = DIO.write_expression(s, d / "expression.txt")
                cfg = DIO.write_settings(d / "settings.txt", expr, tf_file)
            else:
                cfg = d / "settings.txt"
            tasks.append({"arm": arm, "donor": donor, "lineage": lineage,
                          "family": fam, "dir": str(d), "settings": str(cfg),
                          "n_genes": int(s.shape[0]), "n_days": int(s.shape[1]),
                          "days": ",".join(s.columns)})

    rank_rows, run_rows = [], []
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        for res in ex.map(_one_run, tasks):
            rk = res.pop("_ranking")
            run_rows.append(res)
            if rk is not None and not rk.empty:
                rank_rows.append(rk.assign(arm=res["arm"], donor=res["donor"],
                                           lineage=res["lineage"], family=res["family"]))
            log(f"{res['arm']}/{res['donor']}/{res['lineage']}/{res['family']}: "
                f"ok={res['ok']} tfs={res['n_tfs']} ({res['seconds']}s)")

    pd.DataFrame(run_rows).to_csv(OUT / "drem_runs.csv", index=False)
    if rank_rows:
        allrk = pd.concat(rank_rows, ignore_index=True)
        allrk.to_csv(OUT / "drem_tf_rankings.csv", index=False)
        log(f"rankings: {len(allrk)} rows over {allrk.tf.nunique()} TFs")
    json.dump({"n_series": len(series), "n_runs": len(run_rows),
               "n_ok": int(sum(r["ok"] for r in run_rows))},
              open(OUT / "drem_run_summary.json", "w"), indent=1)
    log("done")


if __name__ == "__main__":
    main()
