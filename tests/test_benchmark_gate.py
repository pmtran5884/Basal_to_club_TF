"""Tests for the benchmark gate.

The properties under test are the ones that make the gate meaningful: rare
lineages cannot influence the sweep objective, a control must pass in both
network families, and a TF that ranks equally everywhere fails specificity.
"""
import pandas as pd

from basal_to_club.benchmark.score import (
    check_specificity,
    gating_recovery_score,
    rank_tf,
    score_contrast,
    score_gating,
    score_reported_only,
)

BENCH = {
    "gating": {
        "contrasts": {
            "basal_to_goblet": {"primary": "SPDEF", "secondary": ["FOXA3"],
                                "max_rank_percentile": 0.05, "max_absolute_rank": 25},
            "basal_to_ciliated": {"primary": "FOXJ1", "secondary": ["RFX3"],
                                  "max_rank_percentile": 0.05, "max_absolute_rank": 25},
        },
        "requirements": {"max_fdr": 0.05, "require_positive_nes": True,
                         "require_both_network_families": True,
                         "min_secondary_frac_in_top10pct": 0.5,
                         "specificity": {"SPDEF": {"must_rank_higher_in": "basal_to_goblet",
                                                   "than_in": "basal_to_club"}}},
    },
    "reported_only": {"min_cells": 100, "min_donors": 3,
                      "contrasts": {"basal_to_ionocyte": {"primary": "FOXI1"}}},
}


def activity(top_tfs: dict[str, list[str]], n_filler=400, fdr=1e-6):
    rows = []
    for contrast, ordered in top_tfs.items():
        filler = [f"TF{i}" for i in range(n_filler)]
        for i, tf in enumerate(ordered + filler):
            rows.append({"tf": tf, "contrast": contrast, "nes": 10.0 - 0.01 * i, "fdr": fdr})
    return pd.DataFrame(rows)


def test_rank_is_on_positive_activity():
    """A strongly repressed TF must rank last, not first."""
    df = pd.DataFrame({"tf": ["A", "B"], "contrast": ["c", "c"], "nes": [-9.0, 1.0],
                       "fdr": [1e-9, 0.2]})
    assert rank_tf(df, "c", "A")["rank"] == 2


def test_gate_passes_when_controls_are_recovered_in_both_families():
    act = activity({"basal_to_goblet": ["SPDEF", "FOXA3"],
                    "basal_to_ciliated": ["FOXJ1", "RFX3"],
                    "basal_to_club": ["X"] + [f"F{i}" for i in range(50)] + ["SPDEF"]})
    scores, gate = score_gating({"de_novo": act, "prior": act}, BENCH)
    assert gate["passed"]
    assert gate["specificity_passed"]


def test_gate_fails_when_only_one_family_recovers_the_control():
    good = activity({"basal_to_goblet": ["SPDEF", "FOXA3"],
                     "basal_to_ciliated": ["FOXJ1", "RFX3"],
                     "basal_to_club": ["X", "SPDEF"]})
    bad = activity({"basal_to_goblet": [f"F{i}" for i in range(300)] + ["SPDEF", "FOXA3"],
                    "basal_to_ciliated": ["FOXJ1", "RFX3"],
                    "basal_to_club": ["X", "SPDEF"]})
    _, gate = score_gating({"de_novo": good, "prior": bad}, BENCH)
    assert not gate["passed"]
    assert gate["per_contrast"]["basal_to_goblet"]["families_passed"] == ["de_novo"]


def test_specificity_fails_a_generic_differentiation_signal():
    """SPDEF ranking equally in club and goblet means the signal is not lineage-specific."""
    act = activity({"basal_to_goblet": ["SPDEF"], "basal_to_club": ["SPDEF"]})
    spec = check_specificity(act, BENCH["gating"]["requirements"]["specificity"])
    assert not spec["passed"].all()


def test_rare_lineages_cannot_reach_the_sweep_objective():
    """The objective is computed from gating scores only - structurally, not by convention."""
    act = activity({"basal_to_goblet": ["SPDEF", "FOXA3"],
                    "basal_to_ciliated": ["FOXJ1", "RFX3"],
                    "basal_to_club": ["X"]})
    scores, gate = score_gating({"de_novo": act, "prior": act}, BENCH)
    rare_present = score_reported_only(act, BENCH, {"ionocyte": {"n_cells": 5000, "n_donors": 9}})
    rare_absent = score_reported_only(act, BENCH, {"ionocyte": {"n_cells": 3, "n_donors": 1}})
    assert bool(rare_present.iloc[0]["gates"]) is False
    assert rare_absent.iloc[0]["status"] == "underpowered"
    assert gating_recovery_score(scores) == gate["gating_recovery_score"]


def test_underpowered_rare_lineage_is_recorded_not_omitted():
    act = activity({"basal_to_ionocyte": ["FOXI1"]})
    out = score_reported_only(act, BENCH, {"ionocyte": {"n_cells": 12, "n_donors": 1}})
    assert len(out) == 1
    assert "insufficient" in out.iloc[0]["interpretation"]


def test_missing_tf_fails_rather_than_scoring_nan():
    act = activity({"basal_to_goblet": ["NOTSPDEF"], "basal_to_ciliated": ["FOXJ1", "RFX3"]})
    res = score_contrast(act, "basal_to_goblet", BENCH["gating"]["contrasts"]["basal_to_goblet"],
                         BENCH["gating"]["requirements"])
    assert res[0].passed is False
    assert "absent" in res[0].reason
