"""Signatures with different gene universes must not silently yield all-NaN NES.

aREA propagates a single NaN across a whole sample. Because contrasts are
filtered independently (filterByExpr per contrast), the merged signature matrix
always has NaNs, so passing it whole returned an all-NaN activity matrix with no
error raised. These tests pin the per-signature behaviour with a stub scorer, so
they run without pyviper installed.
"""
import numpy as np
import pandas as pd
import pytest

from basal_to_club.viper.run_viper import apply_per_signature

TFS = ["SPDEF", "FOXJ1", "TP63"]


def _sigs():
    genes = [f"G{i}" for i in range(10)]
    df = pd.DataFrame({"basal_to_club": np.arange(10.0),
                       "basal_to_goblet": np.arange(10.0)[::-1]}, index=genes)
    df.loc["G0", "basal_to_club"] = np.nan       # measured only in one contrast
    df.loc["G9", "basal_to_goblet"] = np.nan
    return df


def _stub(sig: pd.DataFrame) -> pd.DataFrame:
    assert not sig.isna().any().any(), "scorer received NaN genes"
    return pd.DataFrame({sig.columns[0]: np.linspace(1, 3, len(TFS))}, index=TFS)


def test_each_signature_scored_on_its_own_genes():
    out = apply_per_signature(_sigs(), _stub)
    assert list(out.columns) == ["basal_to_club", "basal_to_goblet"]
    assert out.notna().all().all()


def test_all_nan_column_is_an_error_not_a_result():
    def nan_scorer(sig):
        return pd.DataFrame({sig.columns[0]: [np.nan] * len(TFS)}, index=TFS)
    with pytest.raises(ValueError, match="no scores"):
        apply_per_signature(_sigs(), nan_scorer)


def test_empty_signature_is_an_error():
    df = _sigs()
    df["basal_to_ciliated"] = np.nan
    with pytest.raises(ValueError, match="no non-missing genes"):
        apply_per_signature(df, _stub)
