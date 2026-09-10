"""Gene identifier harmonization.

Silent alias mismatches are a real failure mode here: a marker panel that lists
NKX2-1 against a matrix using NKX2_1, or a regulon keyed on an outdated symbol,
degrades quietly rather than erroring. Everything is resolved to current HGNC
symbols once, at load, and the failures are counted.
"""
from __future__ import annotations

import re

import pandas as pd

_SEP = re.compile(r"[._]")


def normalize_symbol(sym: str) -> str:
    """Canonical form used for matching only (NOT written back to var_names)."""
    return _SEP.sub("-", str(sym).strip().upper())


def build_alias_map(hgnc_table: pd.DataFrame) -> dict[str, str]:
    """Map alias/previous symbols -> approved symbol from an HGNC complete set."""
    amap: dict[str, str] = {}
    for _, row in hgnc_table.iterrows():
        approved = row["symbol"]
        amap[normalize_symbol(approved)] = approved
        for col in ("alias_symbol", "prev_symbol"):
            val = row.get(col)
            if isinstance(val, str) and val:
                for alias in val.split("|"):
                    amap.setdefault(normalize_symbol(alias), approved)
    return amap


def resolve(symbols, alias_map: dict[str, str] | None = None) -> list[str]:
    alias_map = alias_map or {}
    return [alias_map.get(normalize_symbol(s), str(s)) for s in symbols]


def match_panel(panel: list[str], available, alias_map: dict[str, str] | None = None):
    """Intersect a marker panel with a matrix's genes. Returns (found, missing).

    Callers must act on `missing`: a club panel that lost SCGB1A1 to an alias
    mismatch will still 'work' and still be wrong.
    """
    alias_map = alias_map or {}
    lookup = {normalize_symbol(g): g for g in available}
    found, missing = [], []
    for gene in panel:
        key = normalize_symbol(alias_map.get(normalize_symbol(gene), gene))
        (found.append(lookup[key]) if key in lookup else missing.append(gene))
    return found, missing


def require_panel(name: str, panel: list[str], available, min_frac: float = 0.5,
                  alias_map: dict[str, str] | None = None) -> list[str]:
    found, missing = match_panel(panel, available, alias_map)
    if len(found) < max(1, int(min_frac * len(panel))):
        raise ValueError(
            f"marker panel '{name}': only {len(found)}/{len(panel)} genes present "
            f"(missing: {missing}). Refusing to score a panel this depleted."
        )
    return found
