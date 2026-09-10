"""Behavioural tests for the club/goblet resolver.

These test the properties the analysis depends on, not the implementation:
goblet cells never become club, gland cells are removed before the call is made,
and a unimodal score distribution does not get a confident threshold.
"""
import numpy as np
import pandas as pd
import pytest

from basal_to_club.secretory.resolver import (
    ResolverThresholds,
    balanced_accuracy,
    classify,
    fit_gmm_threshold,
    fit_thresholds,
    monotonic_club_fraction,
)

RNG = np.random.default_rng(0)


def make_scores(n=400):
    """Two well-separated secretory populations plus contaminants."""
    club = pd.DataFrame({"club": RNG.normal(1.0, 0.2, n), "goblet": RNG.normal(-0.5, 0.2, n),
                         "ciliated": RNG.normal(-0.5, 0.2, n), "smg": RNG.normal(-0.5, 0.2, n)})
    gob = pd.DataFrame({"club": RNG.normal(-0.5, 0.2, n), "goblet": RNG.normal(1.0, 0.2, n),
                        "ciliated": RNG.normal(-0.5, 0.2, n), "smg": RNG.normal(-0.5, 0.2, n)})
    df = pd.concat([club, gob], ignore_index=True)
    df.index = [f"club_{i}" for i in range(n)] + [f"goblet_{i}" for i in range(n)]
    return df


def test_bimodal_threshold_lands_between_modes():
    x = np.concatenate([RNG.normal(0, 0.2, 500), RNG.normal(2, 0.2, 500)])
    thr, sep = fit_gmm_threshold(x)
    assert 0.5 < thr < 1.5
    assert sep > 2


def test_unimodal_scores_do_not_get_a_confident_split():
    """A single Gaussian has no real boundary; separation must report exactly 0.

    This is the case a separation-only criterion gets wrong: a 2-component
    mixture fitted to one Gaussian reports a mean gap around 1.5-2.5 pooled SD.
    """
    thr, sep = fit_gmm_threshold(RNG.normal(0, 1, 2000))
    assert sep == 0.0
    assert thr == pytest.approx(np.quantile(RNG.normal(0, 1, 10), 0.75), abs=10)


def test_club_and_goblet_are_separated():
    scores = make_scores()
    thr = fit_thresholds(scores)
    calls = classify(scores, thr)
    club_correct = (calls[calls.index.str.startswith("club_")] == "club").mean()
    goblet_correct = (calls[calls.index.str.startswith("goblet_")] == "goblet").mean()
    assert club_correct > 0.9
    assert goblet_correct > 0.9


def test_goblet_vetoes_club_when_hybrid_disabled():
    """The asymmetry that protects the club definition from goblet contamination."""
    scores = pd.DataFrame({"club": [2.0], "goblet": [2.0], "ciliated": [-1.0], "smg": [-1.0]})
    thr = ResolverThresholds(club=0.5, goblet=0.5, smg=0.5, ciliated=0.5)
    assert classify(scores, thr, hybrid_class=False, goblet_vetoes_club=True).iloc[0] == "goblet"
    assert classify(scores, thr, hybrid_class=True).iloc[0] == "hybrid"


def test_smg_and_ciliated_override_the_secretory_call():
    """A gland or ciliated cell must never be counted as club."""
    scores = pd.DataFrame({"club": [2.0, 2.0], "goblet": [-1.0, -1.0],
                           "ciliated": [-1.0, 2.0], "smg": [2.0, -1.0]})
    thr = ResolverThresholds(club=0.5, goblet=0.5, smg=0.5, ciliated=0.5)
    calls = classify(scores, thr)
    assert set(calls) == {"smg", "ciliated_contaminant"}


def test_balanced_accuracy_ignores_classes_outside_the_test():
    truth = pd.Series(["club", "goblet", "hybrid", "club"])
    pred = pd.Series(["club", "goblet", "club", "club"])
    assert balanced_accuracy(truth, pred) == pytest.approx(1.0)


def test_club_fraction_rises_over_the_ali_time_course():
    obs = pd.DataFrame({
        "ali_day": np.repeat([0, 7, 14, 21, 28], 100),
        "cell_class": (["basal"] * 100 + ["club"] * 10 + ["basal"] * 90
                       + ["club"] * 30 + ["basal"] * 70 + ["club"] * 60 + ["basal"] * 40
                       + ["club"] * 80 + ["basal"] * 20)})
    out = monotonic_club_fraction(obs)
    assert out["status"] == "ok"
    assert out["monotonic_increase"]


def test_timecourse_check_reports_when_it_cannot_run():
    obs = pd.DataFrame({"ali_day": [np.nan] * 10, "cell_class": ["club"] * 10})
    assert monotonic_club_fraction(obs)["status"] == "no_timecourse_samples"
