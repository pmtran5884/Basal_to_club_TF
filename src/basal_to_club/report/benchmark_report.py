"""Benchmark report. Written whether the gate passed or failed.

A failed gate produces a report explaining WHY it failed - that is a result about
the configuration and is more useful than a missing file.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from basal_to_club.utils.io import setup_logging


def render(gate: dict, scores: pd.DataFrame, netqc: pd.DataFrame,
           resolver: dict) -> str:
    verdict = "PASSED" if gate["passed"] else "FAILED"
    lines = [
        "# Positive-control benchmark", "",
        f"**Gate: {verdict}** — gating recovery score "
        f"`{gate['gating_recovery_score']:.4f}`  ",
        f"Configuration frozen before scoring: `{gate.get('frozen_config')}`", "",
        "The gating controls decide whether the club prediction is produced at all. "
        "The rare-lineage controls below are scored and reported but can never gate, "
        "and never entered the sweep objective.", "",
        "## Gating controls", "",
    ]
    prim = scores[scores["kind"] == "primary"]
    lines += ["| contrast | TF | family | rank | percentile | NES | FDR | pass |",
              "|---|---|---|---|---|---|---|---|"]
    for _, r in prim.iterrows():
        pct = "—" if pd.isna(r["rank_percentile"]) else f"{100*r['rank_percentile']:.2f}%"
        lines.append(f"| {r['contrast']} | **{r['tf']}** | {r.get('family','—')} | "
                     f"{r['rank']} | {pct} | {r['nes']:.2f} | {r['fdr']:.3g} | "
                     f"{'yes' if r['passed'] else 'NO'} |")

    lines += ["", "## Specificity", "",
              "Each control TF must rank higher in its own contrast than in the club "
              "contrast. Without this, a pipeline detecting generic differentiation "
              "would pass.", ""]
    for row in gate.get("specificity_detail", []):
        lines.append(f"- `{row['tf']}` — own contrast rank {row['own_rank']}, "
                     f"club rank {row['club_rank']}: {row['reason']}")

    lines += ["", "## Rare lineages (reported only, never gating)", "",
              "| contrast | TF | status | cells | donors | rank | interpretation |",
              "|---|---|---|---|---|---|---|"]
    for row in gate.get("reported_only", []):
        lines.append(f"| {row['contrast']} | {row['tf']} | {row['status']} | "
                     f"{row['n_cells']} | {row['n_donors']} | {row['rank'] or '—'} | "
                     f"{row['interpretation']} |")

    lines += ["", "## Upstream checks", "",
              f"- Resolver holdout balanced accuracy: "
              f"`{resolver.get('balanced_accuracy')}` "
              f"(threshold {resolver.get('threshold')}, "
              f"passed: {resolver.get('holdout_passed')})",
              f"- ALI time-course club-fraction trend: "
              f"`{resolver.get('ali_timecourse', {}).get('spearman_rho')}`",
              f"- Networks failing regulon sanity: "
              f"`{sorted(netqc.loc[~netqc['passed'], 'network'].unique().tolist())}`", ""]

    if not gate["passed"]:
        lines += ["## Why it failed", ""]
        for contrast, detail in gate["per_contrast"].items():
            if not detail["passed"]:
                lines.append(f"- **{contrast}**: {'; '.join(detail['failure_reasons']) or 'see table'}")
        lines += ["", "The club prediction was not produced. Do not relax "
                  "`config/benchmark.yaml` to clear this.", ""]
    return "\n".join(lines)


def main(sm):
    log = setup_logging(sm.log[0])
    gate = json.load(open(sm.input.gate))
    scores = pd.read_csv(sm.input.scores)
    netqc = pd.read_csv(sm.input.netqc)
    resolver = json.load(open(sm.input.resolver))
    Path(sm.output[0]).write_text(render(gate, scores, netqc, resolver))
    log.info("benchmark report written (gate %s)", "passed" if gate["passed"] else "failed")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
