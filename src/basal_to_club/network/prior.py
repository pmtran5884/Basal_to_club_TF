"""Convert unsigned prior gene-set libraries into signed VIPER regulons.

ChEA3 libraries are sets of TF -> target genes with no direction and no weight.
VIPER needs both: a mode of regulation in [-1, 1] and a likelihood in [0, 1].
Inventing them would make the prior networks unfalsifiable, so both are estimated
from the data being analyzed and the estimate is auditable:

  mode       = sign(rho), where rho is the TF-target Spearman correlation across
               airway metacells (metacells, not single cells, because dropout
               drives single-cell TF-target correlation toward zero)
  likelihood = |rho| scaled by cross-library support - an edge seen in several
               independent ChEA3 libraries is weighted above one seen in a single
               library

Edges below `min_abs_rho` are RETAINED with likelihood 0 rather than dropped.
Dropping them would let the prior network quietly become a data-derived network,
and the whole point of keeping a prior family is that it is not.

The heuristic is audited against CollecTRI, which ships curated signs: see
`audit_signs_against_reference`.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REGULON_COLUMNS = ["tf", "target", "mo", "likelihood"]


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Parse a GMT gene-set library. ChEA3 set names are like 'FOXJ1_ENCODE'."""
    sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split("\t") if p.strip()]
            if len(parts) < 3:
                continue
            name, genes = parts[0], parts[2:]
            tf = name.split("_")[0].split(" ")[0].upper()
            sets.setdefault(tf, []).extend(g.upper() for g in genes)
    return {tf: sorted(set(g)) for tf, g in sets.items()}


def spearman_signs(expr: pd.DataFrame, pairs: list[tuple[str, str]]) -> pd.Series:
    """Spearman rho for each (tf, target) pair over metacells x genes."""
    ranks = expr.rank(axis=0)
    ranks = (ranks - ranks.mean()) / ranks.std(ddof=0).replace(0, np.nan)
    out = {}
    for tf, tgt in pairs:
        if tf in ranks.columns and tgt in ranks.columns:
            a, b = ranks[tf], ranks[tgt]
            rho = float(np.nanmean(a * b))
            out[(tf, tgt)] = 0.0 if not np.isfinite(rho) else float(np.clip(rho, -1, 1))
        else:
            out[(tf, tgt)] = np.nan
    return pd.Series(out, dtype=float)


def build_regulon(gene_sets: dict[str, list[str]], expr: pd.DataFrame,
                  library_support: dict[tuple[str, str], int] | None = None,
                  min_abs_rho: float = 0.10, min_targets: int = 25,
                  max_targets: int = 500) -> pd.DataFrame:
    """Unsigned sets + metacell expression -> signed, weighted regulon table."""
    pairs = [(tf, t) for tf, targets in gene_sets.items() for t in targets
             if t in expr.columns and tf in expr.columns and t != tf]
    if not pairs:
        return pd.DataFrame(columns=REGULON_COLUMNS)

    rho = spearman_signs(expr, pairs).fillna(0.0)
    df = pd.DataFrame({"tf": [p[0] for p in pairs], "target": [p[1] for p in pairs],
                       "rho": rho.to_numpy()})

    df["mo"] = np.sign(df["rho"])
    df.loc[df["mo"] == 0, "mo"] = 1.0          # unresolvable sign -> assume activation
    weak = df["rho"].abs() < min_abs_rho
    df["likelihood"] = df["rho"].abs()
    df.loc[weak, "likelihood"] = 0.0            # kept, but contributes nothing

    if library_support:
        support = df.set_index(["tf", "target"]).index.map(
            lambda k: library_support.get(k, 1))
        df["likelihood"] = df["likelihood"] * np.log1p(np.asarray(support, dtype=float)) / np.log(2)

    # Normalize within TF so regulon size does not become a proxy for confidence.
    df["likelihood"] = df.groupby("tf")["likelihood"].transform(
        lambda s: s / s.max() if s.max() > 0 else s)

    sizes = df.groupby("tf")["target"].size()
    keep = sizes[sizes >= min_targets].index
    df = df[df["tf"].isin(keep)]
    df = (df.sort_values(["tf", "likelihood"], ascending=[True, False])
            .groupby("tf", group_keys=False).head(max_targets))
    return df[REGULON_COLUMNS].reset_index(drop=True)


def cross_library_support(libraries: dict[str, dict[str, list[str]]]) -> dict[tuple[str, str], int]:
    """Count how many independent libraries contain each (tf, target) edge."""
    counts: dict[tuple[str, str], int] = {}
    for sets in libraries.values():
        for tf, targets in sets.items():
            for t in set(targets):
                counts[(tf, t)] = counts.get((tf, t), 0) + 1
    return counts


def audit_signs_against_reference(regulon: pd.DataFrame, reference: pd.DataFrame) -> dict:
    """Compare inferred signs with a curated signed network (CollecTRI).

    Reported, not corrected. If concordance on shared edges is near chance, the
    sign heuristic is not working on this dataset and the prior family should be
    read as unsigned - that is a finding about the build, and it belongs in the
    report rather than in a silent fix.
    """
    merged = regulon.merge(reference[["tf", "target", "mo"]], on=["tf", "target"],
                           suffixes=("", "_ref"))
    if merged.empty:
        return {"n_shared_edges": 0, "concordance": None}
    agree = float((np.sign(merged["mo"]) == np.sign(merged["mo_ref"])).mean())
    weighted = merged[merged["likelihood"] > 0]
    return {
        "n_shared_edges": int(len(merged)),
        "concordance": agree,
        "concordance_weighted_edges_only": (
            float((np.sign(weighted["mo"]) == np.sign(weighted["mo_ref"])).mean())
            if len(weighted) else None),
        "interpretation": ("sign heuristic behaving" if agree >= 0.65 else
                           "near chance - treat prior family as unsigned and say so"),
    }


def regulon_sanity_check(regulon: pd.DataFrame, rules: dict) -> pd.DataFrame:
    """Cheap gate: does FOXJ1 regulate cilia genes, does SPDEF regulate mucins?

    A network that fails this is broken in a way no downstream statistic will
    make obvious, so it is caught before VIPER rather than after.
    """
    rows = []
    for tf, rule in rules.items():
        targets = set(regulon.loc[regulon["tf"] == tf, "target"])
        hits = sorted(targets & set(rule["must_contain_any"]))
        rows.append({"tf": tf, "n_targets": len(targets), "n_hits": len(hits),
                     "hits": ",".join(hits), "min_hits": rule.get("min_hits", 1),
                     "passed": len(hits) >= rule.get("min_hits", 1)})
    return pd.DataFrame(rows)
