"""Configuration invariants.

These are the guardrails that make the design claim true: benchmark thresholds
are frozen, the rare-lineage set never overlaps the gating set, and every
dataset in the registry carries the provenance fields the report needs.
"""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    return yaml.safe_load((ROOT / "config" / f"{name}.yaml").read_text())


def test_gating_and_rare_control_sets_are_disjoint():
    b = load("benchmark")
    gating = set(b["gating"]["contrasts"])
    rare = set(b["reported_only"]["contrasts"])
    assert not (gating & rare), f"a contrast cannot both gate and be reported-only: {gating & rare}"


def test_club_contrast_is_never_a_control():
    b = load("benchmark")
    assert "basal_to_club" not in b["gating"]["contrasts"]
    assert "basal_to_club" not in b["reported_only"]["contrasts"]


def test_every_gating_control_has_a_primary_tf_and_thresholds():
    for name, spec in load("benchmark")["gating"]["contrasts"].items():
        assert spec.get("primary"), name
        assert 0 < spec["max_rank_percentile"] <= 0.2, name
        assert spec["max_absolute_rank"] >= 1, name


def test_datasets_carry_provenance():
    for d in load("datasets")["datasets"]:
        assert d.get("id") and d.get("source") and d.get("role")
        assert d.get("accession") or d.get("cellxgene_dataset_id"), d["id"]


def test_marker_panels_exist_for_every_resolver_program():
    m = load("markers")["markers"]
    for panel in ("club", "goblet", "ciliated", "smg_serous", "smg_mucous", "basal"):
        assert len(m[panel]) >= 3, panel


def test_muc5b_is_not_a_discriminant_marker():
    """MUC5B is expressed by club, goblet and gland mucous cells in human airway."""
    m = load("markers")["markers"]
    assert "MUC5B" not in m["goblet"] or "MUC5B" in m["club"], (
        "MUC5B must not separate club from goblet")
