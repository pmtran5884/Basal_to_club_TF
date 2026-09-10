"""Tests for the unsigned-library to signed-regulon conversion."""
import numpy as np
import pandas as pd

from basal_to_club.network.prior import (
    audit_signs_against_reference,
    build_regulon,
    cross_library_support,
    regulon_sanity_check,
    spearman_signs,
)

RNG = np.random.default_rng(1)


def toy_expression(n=200):
    tf = RNG.normal(size=n)
    return pd.DataFrame({
        "SPDEF": tf,
        "MUC5AC": tf * 0.9 + RNG.normal(0, 0.2, n),      # activated
        "TP63": -tf * 0.9 + RNG.normal(0, 0.2, n),       # repressed
        "NOISE": RNG.normal(size=n),                     # unrelated
        **{f"F{i}": RNG.normal(size=n) for i in range(30)},
    })


def test_sign_follows_the_correlation():
    expr = toy_expression()
    rho = spearman_signs(expr, [("SPDEF", "MUC5AC"), ("SPDEF", "TP63")])
    assert rho[("SPDEF", "MUC5AC")] > 0.5
    assert rho[("SPDEF", "TP63")] < -0.5


def test_weak_edges_are_kept_with_zero_likelihood():
    """Dropping them would turn the prior network into a data-derived one."""
    expr = toy_expression()
    sets = {"SPDEF": ["MUC5AC", "TP63", "NOISE"] + [f"F{i}" for i in range(30)]}
    reg = build_regulon(sets, expr, min_targets=5)
    assert "NOISE" in set(reg["target"])
    assert float(reg.loc[reg["target"] == "NOISE", "likelihood"].iloc[0]) == 0.0
    assert float(reg.loc[reg["target"] == "MUC5AC", "likelihood"].iloc[0]) > 0


def test_repressive_edges_get_negative_mode():
    expr = toy_expression()
    sets = {"SPDEF": ["MUC5AC", "TP63"] + [f"F{i}" for i in range(30)]}
    reg = build_regulon(sets, expr, min_targets=5)
    assert float(reg.loc[reg["target"] == "TP63", "mo"].iloc[0]) == -1.0
    assert float(reg.loc[reg["target"] == "MUC5AC", "mo"].iloc[0]) == 1.0


def test_small_regulons_are_dropped():
    expr = toy_expression()
    reg = build_regulon({"SPDEF": ["MUC5AC", "TP63"]}, expr, min_targets=25)
    assert reg.empty


def test_cross_library_support_counts_independent_evidence():
    libs = {"a": {"SPDEF": ["MUC5AC"]}, "b": {"SPDEF": ["MUC5AC", "TP63"]}}
    support = cross_library_support(libs)
    assert support[("SPDEF", "MUC5AC")] == 2
    assert support[("SPDEF", "TP63")] == 1


def test_sanity_check_flags_a_network_missing_its_canonical_targets():
    reg = pd.DataFrame({"tf": ["FOXJ1"] * 3, "target": ["A", "B", "C"],
                        "mo": 1.0, "likelihood": 1.0})
    rules = {"FOXJ1": {"must_contain_any": ["DNAI1", "TPPP3"], "min_hits": 1}}
    assert not regulon_sanity_check(reg, rules)["passed"].all()


def test_sign_audit_reports_concordance_without_correcting():
    reg = pd.DataFrame({"tf": ["A", "A"], "target": ["X", "Y"], "mo": [1.0, -1.0],
                        "likelihood": [1.0, 1.0]})
    ref = pd.DataFrame({"tf": ["A", "A"], "target": ["X", "Y"], "mo": [1.0, 1.0]})
    audit = audit_signs_against_reference(reg, ref)
    assert audit["n_shared_edges"] == 2
    assert audit["concordance"] == 0.5
