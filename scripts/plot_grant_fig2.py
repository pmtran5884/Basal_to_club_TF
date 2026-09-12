"""Grant figure 2: controls pass, but the two inference methods do not agree on club."""
from pathlib import Path
import json
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "docs/grant_figures"; T = G / "tables"
BLUE, RED, ORANGE, GREEN, GREY, PURPLE = "#4C72B0", "#C44E52", "#DD8452", "#55A868", "#7F7F7F", "#8172B3"

ctrl = pd.read_csv(T / "control_recovery_by_method.csv")
M = pd.read_csv(T / "method_rank_agreement_spearman.csv", index_col=0)
TT = pd.read_csv(T / "method_rank_agreement_top20.csv", index_col=0)
conv = {k: json.load(open(ROOT / f"results/{v}/drem_viper_convergence.json"))
        for k, v in [("nasal", "drem_ali"), ("bronchial", "drem_bronchial")]}
rank_top = pd.read_csv(T / "club_top10_by_method.csv")

fig = plt.figure(figsize=(11.8, 9.0))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.34, left=0.10, right=0.975, top=0.93, bottom=0.08,
                      width_ratios=[1.15, 1.0])

# ---- a: positive-control recovery ------------------------------------------
ax = fig.add_subplot(gs[0, 0])
order = [("ciliated", t) for t in ["FOXJ1", "RFX2", "RFX3", "TP73", "MYB"]] + \
        [("goblet", t) for t in ["SPDEF", "CREB3L1", "XBP1"]] + \
        [("club (prior)", t) for t in ["ELF3", "KLF5", "FOXA2", "CEBPB", "SOX2"]]
ypos = {k: len(order) - 1 - i for i, k in enumerate(order)}
style = {"VIPER static (matched contrast)": (BLUE, "o", 34),
         "DREM nasal (matched lineage)": (RED, "^", 30),
         "DREM bronchial (matched lineage)": (ORANGE, "s", 26)}
for meth, (c, mk, s) in style.items():
    d = ctrl[ctrl.method == meth].dropna(subset=["pct"])
    ax.scatter(d.pct.clip(lower=8e-4), [ypos[(l, t)] for l, t in zip(d.lineage, d.tf)],
               color=c, marker=mk, s=s, label=meth.replace(" (matched contrast)", "").replace(" (matched lineage)", ""), zorder=3)
ax.axvspan(8e-4, 0.05, color=GREEN, alpha=0.08, zorder=0)
ax.axvline(0.05, color=GREEN, lw=0.8, zorder=1)
ax.set_yticks(list(ypos.values()))
short = {"ciliated": "cil", "goblet": "gob", "club (prior)": "club"}
ax.set_yticklabels([f"{t}  ({short[l]})" for l, t in ypos], fontsize=7.5)
for yb in (7.5, 4.5):
    ax.axhline(yb, color=GREY, lw=0.5, ls="-", alpha=0.5, zorder=1)
ax.set_xscale("log"); ax.set_xlim(6e-4, 1.3)
ax.set_xlabel("rank percentile of the control TF in its own lineage contrast\n(left = stronger; shaded = top 5%)")
ax.set_title("a   Controls: ciliated recovered everywhere, goblet only at depth,\n     club priors mid-pack", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="upper right", handletextpad=0.3)

# ---- b: cross-method agreement ---------------------------------------------
ax = fig.add_subplot(gs[0, 1])
labels = [c.replace("\\n", "\n") for c in M.columns]
im = ax.imshow(M.values.astype(float), cmap="RdBu_r", vmin=-1, vmax=1)
for i in range(len(M)):
    for j in range(len(M)):
        v = float(M.values[i, j]); t = int(TT.values[i, j])
        ax.text(j, i, f"{v:.2f}\n{t}/20", ha="center", va="center", fontsize=7,
                color="white" if abs(v) > 0.6 else "black")
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=6.6)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=6.6)
ax.set_title("b   Club rankings: Spearman rho and top-20 overlap between methods", loc="left", fontsize=9)
cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.03); cb.set_label("Spearman rho", fontsize=7.5)
cb.ax.tick_params(labelsize=7)

# ---- c: top-10 lists --------------------------------------------------------
ax = fig.add_subplot(gs[1, 0]); ax.axis("off")
cols = list(rank_top.columns)
colx = np.linspace(0.06, 0.80, len(cols))
for x, c in zip(colx, cols):
    ax.text(x, 1.00, c.replace("\\n", "\n"), fontsize=7.2, ha="left", va="top", fontweight="bold",
            color={0: BLUE, 1: PURPLE, 2: RED, 3: ORANGE}[cols.index(c)], transform=ax.transAxes)
    for i, tf in enumerate(rank_top[c]):
        ax.text(x, 0.86 - i * 0.082, f"{i+1}. {tf}", fontsize=7.2, ha="left", va="top", transform=ax.transAxes)
ax.text(0.0, -0.02, "no TF is shared by a static and a time-course top-10 list", fontsize=7.5,
        color=GREY, transform=ax.transAxes)
ax.set_title("c   Top-10 candidate club drivers, by method", loc="left", fontsize=9)

# ---- d: convergence test ----------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
x = np.arange(2); w = 0.34
vs = [conv[k]["median_pct_viper_set"] for k in ["nasal", "bronchial"]]
bg = [conv[k]["median_pct_background"] for k in ["nasal", "bronchial"]]
ax.bar(x - w / 2, vs, w, color=BLUE, label="VIPER club TF set", zorder=3)
ax.bar(x + w / 2, bg, w, color=GREY, alpha=0.6, label="all other TFs (background)", zorder=3)
ax.axhline(0.5, color="black", lw=0.8, ls="--", zorder=4)
for i, k in enumerate(["nasal", "bronchial"]):
    c = conv[k]
    ax.annotate(f"n={c['n_viper_in_drem']} TFs tested\nMann-Whitney p={c['p']:.2f}\nclub-specific TFs: {c['n_drem_club_specific']}",
                (i, 0.53), ha="center", va="bottom", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(["nasal GSE121600", "bronchial GSE233145"], fontsize=8)
ax.set_ylim(0, 0.78); ax.set_ylabel("median rank percentile in the\nDREM club models (0.5 = chance)")
ax.set_title("d   The time course does not corroborate the static club TF set", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="upper left")

out = G / "grant_fig2_method_disagreement.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
print(out)
