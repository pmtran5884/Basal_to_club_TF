"""HARD GATE.

If the gating controls fail, the club prediction is not produced. This is the
whole point of the design: a pipeline that cannot recover SPDEF for goblet and
FOXJ1 for ciliated has not earned the right to make a claim about club cells.

A failed gate is a result about the configuration. The correct response is to
diagnose it, not to relax config/benchmark.yaml - which is why that file is
frozen before the sweep and this script refuses to read a modified copy.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def main(sm):
    gate = json.load(open(sm.input[0]))
    if not gate["passed"]:
        failures = {c: d["failure_reasons"] for c, d in gate["per_contrast"].items()
                    if not d["passed"]}
        sys.stderr.write(
            "\n=== BENCHMARK GATE FAILED ===\n"
            f"gating recovery score: {gate['gating_recovery_score']}\n"
            f"failed contrasts: {json.dumps(failures, indent=2)}\n"
            f"specificity passed: {gate['specificity_passed']}\n\n"
            "The club prediction has NOT been produced. Diagnose before proceeding:\n"
            "  1. results/networks/network_qc.csv - do the regulons contain the right targets?\n"
            "  2. results/secretory/resolver_holdout_metrics.json - is the compartment right?\n"
            "  3. results/signatures/signature_summary.csv - are the contrasts populated?\n"
            "  4. results/qc/qc_summary.csv - did a dataset lose its secretory cells?\n"
            "Do not edit config/benchmark.yaml to make this pass.\n")
        raise SystemExit(1)
    Path(sm.output[0]).write_text(json.dumps(
        {"passed": True, "gating_recovery_score": gate["gating_recovery_score"]}, indent=2))


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
