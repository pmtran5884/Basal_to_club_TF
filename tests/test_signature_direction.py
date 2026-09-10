"""The limma contrast must be oriented by `case`/`control`, not column order.

This is a regression test for a real inversion: with `~0 + factor(cell_class)`
and a contrast of `colnames[1] - colnames[2]`, R orders the levels
alphabetically, so basal-vs-club produced `basal - club`. Every airway control
(basal vs club, goblet, ciliated) has "basal" sorting first, so all three
signatures came out sign-flipped and SPDEF/FOXJ1 ranked at the BOTTOM.

The test does not need R: it asserts the R source that gets executed encodes
the direction explicitly, which is the property that was violated.
"""
import inspect
import re

from basal_to_club.signature import pseudobulk as PB


def test_contrast_is_not_built_from_column_order():
    src = inspect.getsource(PB.limma_voom_signature)
    assert 'colnames(design)[1]' not in src, (
        "contrast built from alphabetical column order inverts basal-vs-X")


def test_contrast_names_case_and_control_explicitly():
    src = inspect.getsource(PB.limma_voom_signature)
    m = re.search(r'contrasts\s*=\s*"([^"]+)"', src)
    assert m, "expected a literal contrast string"
    lhs, rhs = m.group(1).split("-")
    assert "case" in lhs and "control" in rhs, f"contrast {m.group(1)} is not case-minus-control"


def test_group_column_is_derived_from_case_argument():
    src = inspect.getsource(PB.limma_voom_signature)
    assert "grp_contrast" in src and "== case" in src, "case must define the case level"
