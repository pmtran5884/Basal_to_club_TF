"""Nasal vs bronchial DREM arms: controls, reproducibility, VIPER convergence.

Reads the per-dataset tables written by scripts/compare_drem_viper.py and draws
the four control/comparison panels used in docs/drem_bronchial/.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SETS = {"nasal GSE121600": "drem_ali", "bronchial GSE233145": "drem_bronchial"}
COL = {"nasal GSE121600": "#8C8C8C", "bronchial GSE233145": "#C44E52"}
LINCOL = {"ciliated": "#55A868", "club": "#C44E52", "goblet": "#DD8452"}
GREY = "#8C8C8C"
STYLE = {"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
         "axes.spines.top": False, "axes.spines.right": False,
         "axes.titlelocation": "left", "figure.dpi": 300, "savefig.bbox": "tight"}


def club_pct(setdir: Path) -> pd.DataFrame:
    rk = pd.read_csv(setdir / "drem_tf_rankings.csv")
    rk["pct"] = rk.groupby(["arm", "donor", "lineage", "family"])["rank"].transform(
        lambda r: r / r.max())
    co = rk[(rk.family == "coexpression") & (rk.lineage == "club")]
    p = co.pivot_table(index="tf", columns=["arm", "donor"], values="pct")
    p.columns = [f"{a[:5]}_{d}" for a, d in p.columns]
    # match compare_drem_viper.py: only TFs scored in all four club models
    return p.dropna()


def main():
    ct = {k: pd.read_csv(ROOT / "results" / v / "drem_control_recovery.csv")
          for k, v in SETS.items()}
    rep = {k: pd.read_csv(ROOT / "results" / v / "drem_club_reproducibility.csv")
           for k, v in SETS.items()}
    conv = {k: json.load(open(ROOT / "results" / v / "drem_viper_convergence.json"))
            for k, v in SETS.items()}
    club = {k: club_pct(ROOT / "results" / v) for k, v in SETS.items()}
    vr = pd.read_csv(ROOT / "docs/atlas_label_arms/tables/club_specificity_robustness.csv",
                     index_col=0)
    vset = list(vr[vr.robust_club_specific].index)

    with mpl.rc_context(STYLE):
        fig = plt.figure(figsize=(11, 7.4))
        gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.3)

        # a) control TF recovery, both datasets
        ax = fig.add_subplot(gs[0, 0])
        tabs = {k: v.set_index(["lineage", "tf"])["best"] for k, v in ct.items()}
        idx = sorted(set(tabs["nasal GSE121600"].index) | set(tabs["bronchial GSE233145"].index),
                     key=lambda t: (["ciliated", "club", "goblet"].index(t[0]), t[1]))
        y = np.arange(len(idx))
        for i, key in enumerate(idx):
            a = tabs["nasal GSE121600"].get(key, np.nan)
            b = tabs["bronchial GSE233145"].get(key, np.nan)
            if np.isfinite(a) and np.isfinite(b):
                ax.plot([a, b], [i, i], color=LINCOL[key[0]], lw=0.9, alpha=0.45, zorder=1)
            ax.scatter(a, i, s=24, facecolor="white", edgecolor=LINCOL[key[0]], lw=1.1, zorder=3)
            ax.scatter(b, i, s=26, color=LINCOL[key[0]], zorder=3)
        ax.axvline(0.10, color=GREY, lw=0.8, ls=":")
        ax.annotate("top 10%", (0.115, -0.55), fontsize=6, color=GREY)
        ax.set_yticks(y); ax.set_yticklabels([t for _, t in idx], fontstyle="italic")
        ax.invert_yaxis(); ax.set_xlim(-0.03, 1.03)
        ax.set_title("Lineage control TFs recover in both tissues, except goblet")
        ax.set_xlabel("best rank percentile across the 4 models "
                      "(open = nasal, filled = bronchial; lower = stronger)")
        for lin, c in LINCOL.items():
            ii = [i for i, k in enumerate(idx) if k[0] == lin]
            ax.annotate(lin, (1.0, np.mean(ii)), color=c, fontsize=7, ha="right", va="center")
        ax.text(-0.2, 1.06, "a", transform=ax.transAxes, fontsize=11, fontweight="bold")

        # b) reproducibility
        ax = fig.add_subplot(gs[0, 1])
        w = 0.36
        for i, (k, r) in enumerate(rep.items()):
            ax.bar(np.arange(2) + (i - 0.5) * w, r.mean_spearman, width=w, color=COL[k], label=k)
            for j, v in enumerate(r.mean_spearman):
                ax.annotate(f"{v:.2f}", (j + (i - 0.5) * w, v), xytext=(0, 3),
                            textcoords="offset points", ha="center", fontsize=6.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["same donor,\nboth label arms", "different\ndonors"])
        ax.set_ylabel("mean Spearman of club TF ranks"); ax.set_ylim(0, 1.05)
        ax.set_title("Donor, not the labeling arm, sets the club ranking")
        ax.legend(frameon=False, loc="upper right")
        ax.text(-0.2, 1.06, "b", transform=ax.transAxes, fontsize=11, fontweight="bold")

        # c) VIPER club set vs background in DREM club models
        ax = fig.add_subplot(gs[1, 0])
        pos = 0
        for k, cl in club.items():
            present = [t for t in vset if t in cl.index]
            vv = cl.loc[present].values.ravel(); vv = vv[~np.isnan(vv)]
            bg = cl.drop(index=present).values.ravel(); bg = bg[~np.isnan(bg)]
            U, p = mannwhitneyu(vv, bg, alternative="less")
            assert abs(p - conv[k]["p"]) < 1e-9, (p, conv[k]["p"])
            for j, (v, c) in enumerate([(bg, GREY), (vv, COL[k])]):
                b = ax.violinplot([v], positions=[pos + j], showextrema=False, widths=0.8)
                b["bodies"][0].set_facecolor(c); b["bodies"][0].set_alpha(0.45)
                b["bodies"][0].set_edgecolor("none")
                ax.hlines(np.median(v), pos + j - 0.28, pos + j + 0.28, color="black", lw=1.3, zorder=4)
            ax.annotate(f"p = {p:.2f}", (pos + 0.5, 1.04), ha="center", fontsize=6.5, color=GREY)
            ax.annotate(k.split()[0], (pos + 0.5, 1.13), ha="center", fontsize=7)
            pos += 2.6
        ax.set_xticks([0, 1, 2.6, 3.6])
        ax.set_xticklabels(["all TFs", "VIPER\nclub set"] * 2)
        ax.set_ylim(0, 1.2); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_ylabel("DREM club rank percentile")
        ax.set_title("VIPER's club TFs sit at the median of both time courses")
        ax.text(-0.2, 1.06, "c", transform=ax.transAxes, fontsize=11, fontweight="bold")

        # d) cross-dataset agreement of club rankings
        ax = fig.add_subplot(gs[1, 1])
        n, b = club["nasal GSE121600"], club["bronchial GSE233145"]
        sh = n.index.intersection(b.index)
        xn, xb = n.loc[sh].mean(1), b.loc[sh].mean(1)
        rho = xn.corr(xb, method="spearman")
        ax.scatter(xn, xb, s=9, color=GREY, alpha=0.45, edgecolor="none")
        cons = sorted(set(n.index[(n < 0.10).sum(1) >= 3]) & set(b.index[(b < 0.10).sum(1) >= 3]))
        ax.scatter(xn[cons], xb[cons], s=30, color="#4C72B0", zorder=3)
        off = {"PAX9": (5, 4), "ZBTB7C": (5, -8), "FOXM1": (6, -9), "MYBL2": (6, 5),
               "TFDP1": (-34, 1)}
        for t in cons:
            ax.annotate(t, (xn[t], xb[t]), xytext=off.get(t, (5, 4)),
                        textcoords="offset points", fontsize=6.5, color="#4C72B0",
                        fontstyle="italic")
        for t in [x for x in vset if x in sh]:
            ax.scatter(xn[t], xb[t], s=26, facecolor="none", edgecolor="#C44E52", lw=1.0, zorder=3)
        ax.scatter([], [], s=26, facecolor="none", edgecolor="#C44E52", lw=1.0,
                   label="VIPER robust club set")
        ax.scatter([], [], s=30, color="#4C72B0", label="top decile in both tissues")
        ax.legend(frameon=False, loc="lower right")
        ax.set_xlabel("mean club rank percentile, nasal")
        ax.set_ylabel("mean club rank percentile, bronchial")
        ax.set_title(f"Club rankings agree weakly across tissues ($\\rho$ = {rho:.2f}, n = {len(sh)})")
        ax.text(-0.2, 1.06, "d", transform=ax.transAxes, fontsize=11, fontweight="bold")

        out = ROOT / "results/drem_bronchial/drem_controls_nasal_vs_bronchial.png"
        fig.savefig(out)
        plt.close(fig)
    summary = {"spearman_cross_dataset": round(float(rho), 3), "n_shared_tfs": int(len(sh)),
               "consensus_top_decile_both": cons,
               "top20_overlap": int(len(set(xn.sort_values().head(20).index) &
                                        set(xb.sort_values().head(20).index)))}
    json.dump(summary, open(ROOT / "results/drem_bronchial/cross_dataset_club_agreement.json", "w"),
              indent=1)
    print(json.dumps(summary, indent=1)); print(out)


if __name__ == "__main__":
    main()
