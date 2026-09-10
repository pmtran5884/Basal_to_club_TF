"""Arbitrate three independent opinions about each cell's identity.

Author label, transferred label and marker score are combined by modal vote.
Disagreements are written out rather than resolved silently: where the three
disagree is exactly where a club/goblet mistake would hide.
"""
from __future__ import annotations

import pandas as pd


def marker_call(scores: pd.DataFrame, min_margin: float = 0.05) -> pd.Series:
    """Highest-scoring program per cell, or 'ambiguous' when the top two are close."""
    top2 = scores.apply(lambda r: r.nlargest(2), axis=1)
    best = scores.idxmax(axis=1)
    margin = top2.iloc[:, 0] - top2.iloc[:, 1]
    return best.where(margin >= min_margin, "ambiguous")


def modal_of_three(author: pd.Series, transferred: pd.Series, marker: pd.Series,
                   confidence: pd.Series | None = None,
                   min_confidence: float = 0.6) -> tuple[pd.Series, pd.DataFrame]:
    """Majority vote with an explicit tie-break, plus a disagreement record.

    Tie-break order: transferred (if confident) > marker > author. Author labels
    rank last on ties because they come from heterogeneous annotation vocabularies
    across studies, which is the reason the atlas is being re-annotated at all.
    """
    df = pd.DataFrame({"author": author, "transferred": transferred, "marker": marker})
    if confidence is not None:
        df.loc[confidence < min_confidence, "transferred"] = pd.NA

    def decide(row):
        votes = [v for v in row if pd.notna(v) and v not in ("ambiguous", "unknown",
                                                             "transfer_failed")]
        if not votes:
            return "unassigned"
        counts = pd.Series(votes).value_counts()
        if counts.iloc[0] > 1:
            return counts.index[0]
        for key in ("transferred", "marker", "author"):
            v = row.get(key)
            if pd.notna(v) and v not in ("ambiguous", "unknown", "transfer_failed"):
                return v
        return "unassigned"

    consensus = df.apply(decide, axis=1)
    n_unique = df.apply(lambda r: len({v for v in r if pd.notna(v)}), axis=1)
    disagreement = df.assign(consensus=consensus, n_distinct_opinions=n_unique)
    return consensus, disagreement[disagreement["n_distinct_opinions"] > 1]
