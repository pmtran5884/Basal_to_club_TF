"""Compare a DREM time-course arm against the static VIPER arm.

Runs on either time course (`--drem-set drem_ali` = nasal GSE121600,
`--drem-set drem_bronchial` = bronchial GSE233145) with identical definitions,
so the two datasets are compared on the same footing.

Definitions (fixed here so every dataset is scored the same way):
  pct           within a model (arm x donor x lineage x regulon family), a TF's
                rank divided by the number of TFs scored in that model. Lower is
                stronger; 0.05 means "top 5% of scored TFs".
  control TFs   lineage-specific TFs that must rank high if the method works:
                ciliated FOXJ1/RFX2/RFX3/MYB/TP73, goblet SPDEF/FOXA3/CREB3L1/
                XBP1, club SOX2/FOXA2/CEBPB/NKX2-1/ELF3/KLF5. They are not
                expected to be rank 1, only well inside the head of the list.
  club-specific a TF in the top 10% of the club model and outside the top 25% of
                both the goblet and the ciliated model of the same arm x donor,
                in >=3 of the 4 arm x donor combinations (coexpression family,
                which carries the control TFs; the ChIP-seq family does not
                contain FOXJ1 or SPDEF at all).
  convergence   Mann-Whitney (one-sided) of the DREM club pct of the VIPER
                robust club-specific TFs against all other TFs scored in the
                same club models.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CONTROLS = {"ciliated": ["FOXJ1", "RFX2", "RFX3", "MYB", "TP73"],
            "goblet": ["SPDEF", "FOXA3", "CREB3L1", "XBP1"],
            "club": ["SOX2", "FOXA2", "CEBPB", "NKX2-1", "ELF3", "KLF5"]}
CLUB_TOP = 0.10          # club pct below this
OFF_TOP = 0.25           # goblet and ciliated pct above this
MIN_SUPPORT = 3          # of 4 arm x donor combinations
FAMILY = "coexpression"
STYLE = {"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
         "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
         "axes.spines.top": False, "axes.spines.right": False,
         "axes.titlelocation": "left", "figure.dpi": 300,
         "savefig.bbox": "tight", "font.family": "sans-serif"}
GREY = "#8C8C8C"


def load_rankings(out: Path) -> pd.DataFrame:
    rk = pd.read_csv(out / "drem_tf_rankings.csv")
    rk["pct"] = rk.groupby(["arm", "donor", "lineage", "family"])["rank"].transform(
        lambda r: r / r.max())
    return rk


def control_recovery(rk: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n_tf = rk.groupby(["arm", "donor", "lineage", "family"]).tf.nunique()
    for (a, d, l, f), g in rk.groupby(["arm", "donor", "lineage", "family"]):
        if l == "all":
            continue
        for t in CONTROLS.get(l, []):
            h = g[g.tf == t]
            rows.append({"arm": a, "donor": d, "lineage": l, "family": f, "tf": t,
                         "rank": int(h["rank"].iloc[0]) if len(h) else np.nan,
                         "pct": round(float(h.pct.iloc[0]), 3) if len(h) else np.nan,
                         "n_tf": int(n_tf[(a, d, l, f)])})
    return pd.DataFrame(rows)


def club_specific(piv: pd.DataFrame, combos) -> pd.DataFrame:
    spec = pd.DataFrame(index=piv.index)
    for a, d in combos:
        c, g, ci = [piv.get((a, d, l)) for l in ["club", "goblet", "ciliated"]]
        if c is None or g is None or ci is None:
            continue
        spec[f"{a[:5]}_{d}"] = (c < CLUB_TOP) & (g > OFF_TOP) & (ci > OFF_TOP)
    spec["n_support"] = spec.sum(1)
    spec["pct_club_mean"] = piv[[(a, d, "club") for a, d in combos
                                 if (a, d, "club") in piv.columns]].mean(1).round(3)
    return spec


def reproducibility(piv: pd.DataFrame, combos):
    cl = piv[[(a, d, "club") for a, d in combos if (a, d, "club") in piv.columns]].dropna()
    cl.columns = [f"{a[:5]}_{d}" for a, d, _ in cl.columns]
    rho = cl.corr(method="spearman")
    tops = {c: set(cl[c].sort_values().head(20).index) for c in cl.columns}
    ov = []
    for i, x in enumerate(cl.columns):
        for y in cl.columns[i + 1:]:
            ov.append({"pair": f"{x} vs {y}", "same_donor": x.split("_")[-1] == y.split("_")[-1],
                       "spearman": round(float(rho.loc[x, y]), 3),
                       "overlap_top20": len(tops[x] & tops[y])})
    ov = pd.DataFrame(ov)
    rep = pd.DataFrame({
        "comparison": ["within-donor, across label arm", "across donor"],
        "mean_spearman": [round(float(ov[ov.same_donor].spearman.mean()), 3),
                          round(float(ov[~ov.same_donor].spearman.mean()), 3)],
        "mean_top20_overlap": [round(float(ov[ov.same_donor].overlap_top20.mean()), 1),
                               round(float(ov[~ov.same_donor].overlap_top20.mean()), 1)]})
    return cl, ov, rep


def figure(ct: pd.DataFrame, cl: pd.DataFrame, rep: pd.DataFrame, vrank: pd.DataFrame,
           conv: dict, label: str, path: Path):
    with mpl.rc_context(STYLE):
        fig = plt.figure(figsize=(10.5, 3.6))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.42)
        col = {"ciliated": "#55A868", "goblet": "#DD8452", "club": "#C44E52"}

        ax = fig.add_subplot(gs[0, 0])
        tab = ct[ct.family == FAMILY].pivot_table(index=["lineage", "tf"], values="pct",
                                                  aggfunc=["min", "median"])
        tab.columns = ["best", "median"]
        tab = tab.sort_values(["lineage", "median"])
        y = np.arange(len(tab))
        for i, ((lin, tf), r) in enumerate(tab.iterrows()):
            ax.plot([r["best"], r["median"]], [i, i], color=col[lin], lw=1.1, alpha=0.5)
            ax.scatter(r["median"], i, s=26, color=col[lin], zorder=3)
            ax.scatter(r["best"], i, s=16, facecolor="white", edgecolor=col[lin], lw=1, zorder=3)
        ax.axvline(0.10, color=GREY, lw=0.8, ls=":")
        ax.set_yticks(y); ax.set_yticklabels([f"{tf}" for _, tf in tab.index], fontstyle="italic")
        ax.invert_yaxis(); ax.set_xlim(-0.03, 1.03)
        ax.set_xlabel("rank percentile in its lineage model (lower = stronger)")
        ax.set_title(f"Lineage control TFs, {label}")
        ax.annotate("open = best model\nfilled = median of 4", (0.62, 0.12),
                    xycoords="axes fraction", fontsize=6, color=GREY)
        ax.annotate("top 10%", (0.105, len(tab) - 0.6), fontsize=6, color=GREY)
        for lin, c in col.items():
            ax.plot([], [], color=c, lw=3, label=lin)
        ax.legend(frameon=False, loc="lower right", fontsize=6.5, handlelength=1.2)
        ax.text(-0.28, 1.06, "a", transform=ax.transAxes, fontsize=11, fontweight="bold")

        ax = fig.add_subplot(gs[0, 1])
        x = np.arange(2)
        ax.bar(x, rep.mean_spearman, width=0.55, color=["#4C72B0", GREY])
        for i, v in enumerate(rep.mean_spearman):
            ax.annotate(f"{v:.2f}", (i, v), xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=7)
        ax.set_xticks(x); ax.set_xticklabels(["same donor,\nboth label arms", "different\ndonors"])
        ax.set_ylabel("mean Spearman of club TF ranks")
        ax.set_ylim(0, 1.05); ax.set_title("Donor, not label arm, sets the club ranking")
        ax.text(-0.3, 1.06, "b", transform=ax.transAxes, fontsize=11, fontweight="bold")

        ax = fig.add_subplot(gs[0, 2])
        bg = cl.drop(index=[t for t in vrank.index if t in cl.index]).values.ravel()
        bg = bg[~np.isnan(bg)]
        vv = vrank.values.ravel(); vv = vv[~np.isnan(vv)]
        parts = ax.violinplot([bg, vv], showextrema=False, widths=0.8)
        for b, c in zip(parts["bodies"], [GREY, "#C44E52"]):
            b.set_facecolor(c); b.set_alpha(0.45); b.set_edgecolor("none")
        for i, v in enumerate([bg, vv], start=1):
            ax.hlines(np.median(v), i - 0.28, i + 0.28, color="black", lw=1.4, zorder=4)
        ax.set_xticks([1, 2]); ax.set_xticklabels([f"all TFs\n(n={len(bg)})",
                                                   f"VIPER club set\n(n={len(vv)})"])
        ax.set_ylabel("DREM club rank percentile")
        ax.set_title("VIPER's club TFs are not enriched in DREM")
        ax.annotate(f"one-sided p = {conv['p']:.2f}", (0.5, 0.04), xycoords="axes fraction",
                    ha="center", fontsize=7, color=GREY)
        ax.text(-0.3, 1.06, "c", transform=ax.transAxes, fontsize=11, fontweight="bold")
        fig.savefig(path)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drem-set", default="drem_bronchial")
    ap.add_argument("--viper-robust", default="docs/atlas_label_arms/tables/club_specificity_robustness.csv")
    ap.add_argument("--label", default=None, help="dataset name for figure titles")
    args = ap.parse_args()
    out = ROOT / "results" / args.drem_set
    label = args.label or args.drem_set.replace("drem_", "")

    rk = load_rankings(out)
    ct = control_recovery(rk)
    ctab = ct[ct.family == FAMILY].pivot_table(index=["lineage", "tf"],
                                               columns=["arm", "donor"], values="pct")
    ctab.columns = [f"{a[:5]}_{d}" for a, d in ctab.columns]
    ctab["best"] = ctab.min(1).round(3)
    ctab.round(3).to_csv(out / "drem_control_recovery.csv")
    ct.to_csv(out / "drem_control_recovery_long.csv", index=False)

    co = rk[rk.family == FAMILY]
    piv = co.pivot_table(index="tf", columns=["arm", "donor", "lineage"], values="pct")
    combos = sorted({(a, d) for a, d, _ in piv.columns})
    spec = club_specific(piv, combos)
    spec.sort_values("pct_club_mean").to_csv(out / "drem_club_specific.csv")
    cl, ov, rep = reproducibility(piv, combos)
    rep.to_csv(out / "drem_club_reproducibility.csv", index=False)
    ov.to_csv(out / "drem_club_pairwise.csv", index=False)

    vr = pd.read_csv(ROOT / args.viper_robust, index_col=0)
    vset = [t for t in vr[vr.robust_club_specific].index if t in cl.index]
    vrank = cl.loc[vset]
    bg = cl.drop(index=vset).values.ravel(); bg = bg[~np.isnan(bg)]
    vv = vrank.values.ravel(); vv = vv[~np.isnan(vv)]
    U, p = mannwhitneyu(vv, bg, alternative="less")
    conv = {"drem_set": args.drem_set, "n_viper_robust": int(vr.robust_club_specific.sum()),
            "n_viper_in_drem": len(vset), "median_pct_viper_set": round(float(np.median(vv)), 3),
            "median_pct_background": round(float(np.median(bg)), 3),
            "n_obs_viper": int(len(vv)), "n_obs_background": int(len(bg)),
            "U": float(U), "p": float(p),
            "n_drem_club_specific": int((spec.n_support >= MIN_SUPPORT).sum()),
            "drem_club_specific": sorted(spec[spec.n_support >= MIN_SUPPORT].index),
            "overlap_with_viper": sorted(set(spec[spec.n_support >= MIN_SUPPORT].index) & set(vset))}
    json.dump(conv, open(out / "drem_viper_convergence.json", "w"), indent=1)
    figure(ct, cl, rep, vrank, conv, label, out / f"drem_vs_viper_{label}.png")

    print(ctab.round(3).to_string())
    print("\n", rep.to_string(index=False))
    print("\nconvergence:", json.dumps({k: v for k, v in conv.items()
                                        if k != "drem_club_specific"}, indent=1))


if __name__ == "__main__":
    main()
