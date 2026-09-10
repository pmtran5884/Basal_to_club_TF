"""Final report: the club prediction with its supporting evidence."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from basal_to_club.utils.io import load_yaml, setup_logging


def render(ranking: pd.DataFrame, stability: pd.DataFrame, gate: dict,
           counts: dict, datasets: list[dict]) -> str:
    top = ranking.head(30).join(stability[["median_rank", "frac_runs_in_top_n"]],
                                how="left")
    lines = [
        "# Candidate transcription factors driving basal → club",
        "",
        f"Gate passed with recovery score `{gate['gating_recovery_score']:.4f}`. "
        "The ranking below was produced under the frozen configuration; no "
        "parameter was chosen using a club result.",
        "",
        "## Composition of the resolved atlas",
        "",
        "| cell class | cells | donors |", "|---|---|---|",
    ]
    for cls, c in sorted(counts.items(), key=lambda kv: -kv[1]["n_cells"]):
        lines.append(f"| {cls} | {c['n_cells']} | {c['n_donors']} |")

    lines += ["", "## Ranked candidates", "",
              "`tier` is high_confidence when the TF is positive in both network "
              "families, supported by at least two signature constructions, and "
              "FDR < 0.05. `prior_plausible` is an annotation only — it was not "
              "used in the ranking.", "",
              "| rank | TF | mean NES | de novo | prior | min FDR | methods | tier | "
              "median rank across resamples | prior plausible |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for tf, r in top.iterrows():
        lines.append(
            f"| {int(r['rank'])} | **{tf}** | {r['mean_nes']:.2f} | "
            f"{r.get('de_novo', float('nan')):.2f} | {r.get('prior', float('nan')):.2f} | "
            f"{r['min_fdr']:.3g} | {int(r['n_methods_supporting'])} | {r['tier']} | "
            f"{r.get('median_rank', float('nan'))} | "
            f"{'yes' if r['prior_plausible_club_tf'] else ''} |")

    lines += ["", "## How to read this", "",
              "- A high-confidence TF is a *nomination*, not a demonstrated driver. "
              "VIPER infers activity from the expression of a TF's targets; a TF whose "
              "regulon overlaps the club program will score highly whether or not it is "
              "causal in this system.",
              "- Candidates supported by only one network family are listed but should "
              "be treated as network-specific until corroborated.",
              "- Robustness columns are empty when resampling was not run; in that case "
              "stability is unassessed rather than good.",
              "", "## Datasets used", "", "| id | accession | role |", "|---|---|---|"]
    for d in datasets:
        lines.append(f"| {d['id']} | {d.get('accession','—')} | {d.get('role','—')} |")
    return "\n".join(lines)


def main(sm):
    log = setup_logging(sm.log[0])
    ranking = pd.read_csv(sm.input.ranking, index_col=0)
    stability = pd.read_csv(sm.input.stability, index_col=0)
    gate = json.load(open(sm.input.gate))
    datasets = load_yaml("config/datasets.yaml")["datasets"]
    Path(sm.output[0]).write_text(
        render(ranking, stability, gate, gate.get("lineage_cell_counts", {}), datasets))
    log.info("final report written: %d candidates, %d high-confidence",
             len(ranking), int((ranking["tier"] == "high_confidence").sum()))


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
