"""Grant figure 1: what the public data can and cannot support."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "docs/grant_figures"; T = G / "tables"
BLUE, RED, ORANGE, GREEN, GREY = "#4C72B0", "#C44E52", "#DD8452", "#55A868", "#7F7F7F"

land = pd.read_csv(T / "dataset_landscape.csv")
meth = pd.read_csv(T / "club_counts_by_method.csv")
cov = pd.read_csv(T / "ali_timecourse_coverage.csv")
pb = pd.read_csv(ROOT / "docs/qc_2026-09/run_tables/pseudobulk_meta.csv")

fig = plt.figure(figsize=(11.6, 8.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30, left=0.085, right=0.975, top=0.93, bottom=0.085)

# ---- a: landscape -----------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
ax.axhspan(2.5, 9, color=GREEN, alpha=0.07, zorder=0)
ax.axvspan(3000, 4e4, ymin=0, ymax=1, color=GREEN, alpha=0.07, zorder=0)
atl = land[land.kind.str.startswith("tissue")]
ali = land[~land.kind.str.startswith("tissue")]
ax.scatter(atl.club_published.clip(lower=3), atl.timepoints + np.random.default_rng(1).uniform(-.12, .12, len(atl)),
           s=18 + 6 * atl.donors, facecolor="white", edgecolor=BLUE, lw=1.3, zorder=3, label="tissue atlas (HLCA), no time axis")
ax.scatter(ali.club_pipeline, ali.timepoints, s=18 + 6 * ali.donors, color=RED, zorder=4,
           label="ALI differentiation time course")
for (_, r), dy in zip(ali.iterrows(), (-24, 14)):
    ax.annotate(f"{r.dataset.replace(' ALI ', chr(10))}: {int(r.donors)} donors, {int(r.timepoints)} days",
                (r.club_pipeline, r.timepoints), textcoords="offset points", xytext=(-10, dy),
                ha="right", fontsize=6.8, color=RED)
ax.annotate("10 HLCA tissue studies\n(4-7,095 club cells each,\nno differentiation time axis)",
            (60, 0), textcoords="offset points", xytext=(0, 26), ha="center", fontsize=6.8, color=BLUE)
ax.text(1.1e4, 5.0, "required for a\ntrajectory-based\nTF inference", fontsize=7, color=GREEN, ha="center", va="center")
ax.set_xscale("log"); ax.set_xlim(2, 4e4); ax.set_ylim(-0.9, 9)
ax.set_xlabel("club cells in the dataset (log scale)")
ax.set_ylabel("sampled differentiation time points")
ax.set_title("a   No public dataset has both club depth and a time axis", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="upper left", handletextpad=0.3)

# ---- b: club counts by method ----------------------------------------------
ax = fig.add_subplot(gs[0, 1])
m = meth.sort_values("published_atlas", ascending=True)
y = np.arange(len(m)); h = 0.26
for off, col, c, lab in [(-h, "published_atlas", BLUE, "published atlas annotation"),
                         (0.0, "pipeline_class", ORANGE, "HLCA label transfer (this pipeline)"),
                         (h, "resolver", RED, "marker-based club/goblet resolver")]:
    ax.barh(y + off, m[col].clip(lower=0.7), height=h, color=c, label=lab, zorder=3)
ax.set_yticks(y); ax.set_yticklabels([s.replace("_", " ") for s in m.study], fontsize=6.8)
ax.set_xscale("log"); ax.set_xlim(0.7, 2e4)
ax.set_xlabel("cells called club (log scale)")
ax.set_title("b   The same cells, three club definitions: up to 25-fold disagreement", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="lower right")
for _, r in m.iterrows():
    if r.fold_disagreement >= 5:
        ax.annotate(f"{r.fold_disagreement:.0f}x", (max(r.published_atlas, r.pipeline_class, r.resolver) * 1.5,
                    y[list(m.study).index(r.study)]), fontsize=6.5, va="center", color=GREY)

# ---- c: club cells per donor entering the inference ------------------------
ax = fig.add_subplot(gs[1, 0])
club = pb[pb.cell_class == "club"]["n_cells"].values
cil = pb[pb.cell_class == "ciliated"]["n_cells"].values
bins = np.logspace(np.log10(20), np.log10(4000), 22)
ax.hist(cil, bins=bins, color=GREY, alpha=0.35, label=f"ciliated (n={len(cil)} donor-samples)", zorder=2)
ax.hist(club, bins=bins, color=RED, alpha=0.75, label=f"club (n={len(club)} donor-samples)", zorder=3)
ax.axvline(np.median(club), color=RED, ls="--", lw=1.2, zorder=4)
ax.axvline(25, color="black", lw=1.0, zorder=4)
ax.annotate(f"median {np.median(club):.0f} club cells\nper donor-sample", (np.median(club), ax.get_ylim()[1] * 0.78),
            xytext=(8, 0), textcoords="offset points", fontsize=7, color=RED)
ax.annotate("pseudobulk floor\n(25 cells)", (25, ax.get_ylim()[1] * 0.45), xytext=(6, 0),
            textcoords="offset points", fontsize=6.8)
ax.set_xscale("log"); ax.set_xlabel("cells per donor-sample entering the contrast (log scale)")
ax.set_ylabel("donor-samples")
ax.set_title("c   Club is the scarce arm in every study", loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="upper right")

# ---- d: ALI day coverage ----------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
cov = cov.copy(); cov["row"] = cov.dataset + " | " + cov.donor
rows = sorted(cov.row.unique(), reverse=True)
ypos = {r: i for i, r in enumerate(rows)}
for _, r in cov.iterrows():
    ax.scatter(r.day, ypos[r.row], s=8 + r.club_atlas / 6, color=RED if "bronch" in r.dataset else BLUE,
               alpha=0.85, zorder=3)
    dy = 9 if (r.day // 2) % 2 == 0 else -13
    ax.annotate(f"{int(r.club_atlas)}", (r.day, ypos[r.row]), textcoords="offset points", xytext=(0, dy),
                ha="center", fontsize=5.8, color=GREY)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r.replace(" GSE", "\nGSE").replace(" | ", "  ") for r in rows], fontsize=6.5)
ax.set_xlabel("days at air-liquid interface")
ax.set_title("d   Two donors per time course, and no replicated nasal day grid", loc="left", fontsize=9)
ax.set_xlim(-3, 51); ax.set_ylim(-0.7, len(rows) - 0.3)
ax.annotate("point size and number = club cells at that day", (0.99, 0.03), xycoords="axes fraction",
            ha="right", fontsize=6.5, color=GREY)
for d in [7, 12]:
    ax.plot([d, d], [-0.6, 1.6], color=BLUE, lw=0.7, ls=":", zorder=1)
ax.annotate("only 2 of the 8 nasal days\nare sampled in both donors", (15, 0.55), fontsize=6.5,
            color=BLUE, ha="left", va="center")

out = G / "grant_fig1_data_landscape.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
print(out)
