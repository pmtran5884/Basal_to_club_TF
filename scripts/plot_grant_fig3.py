"""Grant figure 3: what limits the inference, and what a new ALI time course must deliver."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "docs/grant_figures"; T = G / "tables"
BLUE, RED, ORANGE, GREEN, GREY = "#4C72B0", "#C44E52", "#DD8452", "#55A868", "#7F7F7F"

sh = pd.read_csv(T / "power_curve_summary.csv")
lad = pd.read_csv(T / "agreement_ladder.csv")
sb = pd.read_csv(T / "donor_reliability_spearman_brown.csv")
dep = pd.read_csv(ROOT / "results/drem_bronchial/platform_depth_control.csv")

fig = plt.figure(figsize=(11.8, 8.8))
gs = fig.add_gridspec(2, 2, hspace=0.44, wspace=0.30, left=0.09, right=0.975, top=0.93, bottom=0.09)

# ---- a: rank accuracy vs number of case cells ------------------------------
ax = fig.add_subplot(gs[0, 0])
for lin, c, mk in [("club", RED, "o"), ("ciliated", BLUE, "s")]:
    d = sh[sh.lineage == lin].sort_values("n_cells")
    ax.plot(d.n_cells, d.top20_vs_full / 20, color=c, marker=mk, ms=4.5, lw=1.5, label=f"{lin}", zorder=3)
ax.axvspan(40, 100, color=GREY, alpha=0.18, zorder=0)
ax.annotate("contrast cannot\nbe formed at all\n(no donor-day group\nclears 25 cells)", (43, 0.42), fontsize=6.4, color=GREY, va="center")
ax.axhline(0.95, color=GREEN, lw=0.8, ls="--", zorder=2)
ax.annotate("19/20 top TFs recovered", (300, 0.957), fontsize=6.8, color=GREEN, ha="left", va="bottom")
ax.axvline(88, color=RED, lw=1.0, ls=":", zorder=2)
ax.annotate("median public\ndataset: 88 club\ncells per donor", (95, 0.60), fontsize=6.8, color=RED, va="center")
ax.set_xscale("log"); ax.set_xlim(40, 1800); ax.set_ylim(0.25, 1.03)
ax.set_xlabel("club (or ciliated) cells in the case arm of the contrast")
ax.set_ylabel("top-20 TF concordance with the\nfull-coverage ranking")
ax.set_title("a   Rank accuracy versus cells available: club needs ~7x more\n     cells than ciliated to stabilise", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7.5, loc="lower right", title="lineage", title_fontsize=7.5)

# ---- b: agreement ladder ----------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
lab, val, col = [], [], []
for lvl in lad.level.unique():
    d = lad[lad.level == lvl]
    for _, r in d.iterrows():
        lab.append(lvl.replace("\\n", "\n") + ("" if r.dataset == "both" else f"\n[{r.dataset}]"))
        val.append(r.rho)
        col.append(GREY if "technical" in lvl else RED if "labeling" in lvl else ORANGE if "donor" in lvl else BLUE if "tissue" in lvl else "#8172B3")
y = np.arange(len(lab))[::-1]
ax.barh(y, val, color=col, height=0.62, zorder=3)
for yy, v in zip(y, val):
    ax.text(v + 0.015, yy, f"{v:.2f}", va="center", fontsize=7)
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=6.6)
ax.set_xlim(0, 1.08); ax.set_xlabel("Spearman rho between two club TF rankings")
ax.set_title("b   Sampling is not the bottleneck: donor, cell definition and\n     inference method are", loc="left", fontsize=9)

# ---- c: donors required -----------------------------------------------------
ax = fig.add_subplot(gs[1, 0])
for col_, c, lab_ in [("reliability_bronchial_r0.493", BLUE, "observed bronchial inter-donor rho = 0.49"),
                      ("reliability_nasal_r0.228", RED, "observed nasal inter-donor rho = 0.23")]:
    ax.plot(sb.k_donors, sb[col_], color=c, lw=1.8, zorder=3, label=lab_)
ax.axhline(0.8, color=GREEN, lw=0.9, ls="--", zorder=2)
ax.annotate("target reliability 0.8", (24, 0.81), fontsize=6.8, color=GREEN, ha="right")
for col_, c, k in [("reliability_bronchial_r0.493", BLUE, 4.1), ("reliability_nasal_r0.228", RED, 13.5)]:
    ax.plot([k, k], [0, 0.8], color=c, lw=0.9, ls=":", zorder=2)
    ax.annotate(f"{np.ceil(k):.0f} donors", (k, 0.22 if k > 8 else 0.10), color=c, fontsize=7, ha="center",
                bbox=dict(fc="white", ec="none", pad=0.6))
ax.axvspan(0, 2, color=GREY, alpha=0.18, zorder=0)
ax.annotate("every existing\nALI time course:\n2 donors", (2.4, 0.62), fontsize=6.8, color=GREY)
ax.set_xlim(0, 24); ax.set_ylim(0, 1.0)
ax.set_xlabel("donors in the time course (k)")
ax.set_ylabel("reliability of the k-donor consensus TF ranking\n(Spearman-Brown)")
ax.set_title("c   Donors required for a reproducible club TF ranking", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="lower right")

# ---- d: design specification -----------------------------------------------
ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
med_b = int(dep.loc[dep.dataset == "bronchial_dropseq", "median_genes_per_cell"].iloc[0])
med_n = int(dep.loc[dep.dataset == "nasal_10x", "median_genes_per_cell"].iloc[0])
muc_b = float(dep.loc[dep.dataset == "bronchial_dropseq", "pct_MUC5AC_detected"].iloc[0])
rows = [("design parameter", "best available public data", "required (this proposal)"),
        ("donors on one common day grid", "2", "8-12"),
        ("time points sampled in every donor", "nasal 2 of 8; bronchial 7", ">=6"),
        ("club cells per donor-level contrast", "median 88", ">=1,500 (~250 per donor-day)"),
        ("median genes per cell", f"{med_b} bronchial / {med_n:,} nasal", ">=3,000"),
        ("MUC5AC+ cells (goblet exclusion)", f"{muc_b:.2f}% (bronchial)", ">=5%"),
        ("club identity", "3 label definitions disagree\nup to 25-fold", "marker-resolved club vs\ngoblet at adequate depth"),
        ("tissue", "nasal ALI, or bronchial\nat low depth", "bronchial HBEC ALI")]
yy = np.linspace(0.93, 0.06, len(rows))
for i, (a, b, c) in enumerate(rows):
    w = "bold" if i == 0 else "normal"
    fs = 7.2 if i == 0 else 6.7
    ax.text(0.00, yy[i], a, fontsize=fs, fontweight=w, va="center", transform=ax.transAxes)
    ax.text(0.42, yy[i], b, fontsize=fs, fontweight=w, va="center", color=GREY if i else "black", transform=ax.transAxes)
    ax.text(0.745, yy[i], c, fontsize=fs, fontweight=w, va="center", color=GREEN if i else "black", transform=ax.transAxes)
    if i == 0:
        ax.plot([0, 1.0], [yy[i] - 0.045] * 2, color="black", lw=0.8, transform=ax.transAxes, clip_on=False)
ax.set_title("d   Specification for the proposed HBEC ALI time course", loc="left", fontsize=9)

out = G / "grant_fig3_power_and_design.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
print(out)
