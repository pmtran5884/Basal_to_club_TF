"""Parse DREM batch output into a TF ranking.

DREM writes one file per out-edge of every split:
    path_<depth>_<parent_mean>_<parent_sigma>_to_<child_mean>_<child_sigma>.txt
whose rows are TFs with enrichment statistics for the genes taking that edge.
"Score Split" is the p-value against the sibling branches at the same split;
"Score Overall" is against all genes. We rank on the split score where present.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

_PATH_RE = re.compile(
    r"path_(?P<depth>\d+)_(?P<pmean>-?[\d.]+)_(?P<psigma>-?[\d.]+)_to_"
    r"(?P<cmean>-?[\d.]+)_(?P<csigma>-?[\d.]+)\.txt$"
)


def parse_edges(tf_dir: str | Path) -> pd.DataFrame:
    """Long table: one row per (TF, edge) with edge geometry and scores."""
    frames = []
    for p in sorted(Path(tf_dir).glob("path_*.txt")):
        m = _PATH_RE.search(p.name)
        if m is None:
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except pd.errors.EmptyDataError:   # DREM writes a header-less file for
            continue                       # a branch with no scored regulators
        df.columns = [c.strip() for c in df.columns]
        g = m.groupdict()
        df["depth"] = int(g["depth"])
        df["parent_mean"] = float(g["pmean"])
        df["child_mean"] = float(g["cmean"])
        df["delta"] = df["child_mean"] - df["parent_mean"]
        df["edge"] = p.stem
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out = out.rename(columns={out.columns[0]: "tf"})
    # DREM suffixes the TF name with its input level ("FOXM1 1")
    out["tf"] = out.tf.astype(str).str.replace(r"\s+\d+$", "", regex=True)
    return out


def rank_tfs(edges: pd.DataFrame, direction: str = "up",
             min_targets: int = 5) -> pd.DataFrame:
    """Rank TFs by their most significant enrichment on an edge of the given
    direction. `direction='up'` keeps edges whose child mean exceeds the parent,
    i.e. the branch on which the program is induced."""
    if edges.empty:
        return pd.DataFrame(columns=["tf", "score", "rank"])
    df = edges.copy()
    if direction == "up":
        df = df[df.delta > 0]
    elif direction == "down":
        df = df[df.delta < 0]
    score_col = "Score Split" if "Score Split" in df.columns else "Score Overall"
    num_col = "Num Path" if "Num Path" in df.columns else None
    if num_col is not None:
        df = df[pd.to_numeric(df[num_col], errors="coerce").fillna(0) >= min_targets]
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col])
    if df.empty:
        return pd.DataFrame(columns=["tf", "score", "rank"])
    best = (df.sort_values(score_col)
              .groupby("tf", as_index=False)
              .first()
              .rename(columns={score_col: "score"}))
    best["n_edges"] = df.groupby("tf").size().reindex(best.tf).values
    best = best.sort_values("score").reset_index(drop=True)
    best["rank"] = np.arange(1, len(best) + 1)
    keep = ["tf", "score", "rank", "n_edges", "delta", "depth", "edge"]
    return best[[c for c in keep if c in best.columns]]


def read_gene_assignment(path: str | Path) -> pd.DataFrame:
    """geneassign.txt: gene, encoded path id, then the node mean at each step."""
    rows = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        rows.append({"gene": f[0], "path_id": f[1],
                     "means": [float(x) for x in f[2:] if x != ""]})
    df = pd.DataFrame(rows)
    if not df.empty:
        df["final_mean"] = df.means.apply(lambda m: m[-1] if m else np.nan)
    return df


_SPLIT2_RE = re.compile(
    r"split_(?P<depth>\d+)_(?P<pmean>-?[\d.]+)_(?P<psigma>[\d.]+)_2way\.txt$")
_SPLITV_RE = re.compile(
    r"split_(?P<depth>\d+)_(?P<pmean>-?[\d.]+)_(?P<psigma>[\d.]+)_to_"
    r"(?P<cmean>-?[\d.]+)_(?P<csigma>[\d.]+)_vsOthers\.txt$")


def parse_splits(tf_dir: str | Path) -> pd.DataFrame:
    """Long table from DREM's split tables: one row per (TF, branch).

    Two file flavours are written. `_2way` compares the two children of a
    binary split ("Low"/"High" by node mean); `_vsOthers` compares one child of
    a >=3-way split against its siblings. Both give a logistic coefficient, the
    counts of bound/unbound genes either side, the mean TF-input on each side
    and a p-value ("Score").
    """
    frames = []
    for p in sorted(Path(tf_dir).glob("split_*.txt")):
        m2, mv = _SPLIT2_RE.search(p.name), _SPLITV_RE.search(p.name)
        m = m2 or mv
        if m is None:
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        df.columns = [c.strip() for c in df.columns]
        g = m.groupdict()
        tgt, oth = ("High", "Low") if m2 else ("High", "Other")
        out = pd.DataFrame({
            "tf": df[df.columns[0]].astype(str).str.replace(r"\s+\d+$", "", regex=True),
            "coeff": pd.to_numeric(df.get("Coeff"), errors="coerce"),
            "n_target_bound": pd.to_numeric(df.get(f"{tgt} 1"), errors="coerce"),
            "avg_target": pd.to_numeric(df.get(f"Avg. {tgt}"), errors="coerce"),
            "avg_other": pd.to_numeric(df.get(f"Avg. {oth}"), errors="coerce"),
            "score": pd.to_numeric(df.get("Score"), errors="coerce"),
        })
        out["depth"] = int(g["depth"])
        out["parent_mean"] = float(g["pmean"])
        out["child_mean"] = float(g["cmean"]) if mv else np.nan
        out["kind"] = "vsOthers" if mv else "2way"
        out["split"] = p.stem
        frames.append(out)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df["delta"] = df.child_mean - df.parent_mean
    return df


def rank_tfs_from_splits(splits: pd.DataFrame, direction: str = "up",
                         min_targets: int = 5) -> pd.DataFrame:
    """Rank TFs by their most significant association with an induced branch.

    A TF qualifies on a branch when its targets are over-represented among the
    genes taking that branch (`avg_target > avg_other`) and the branch is the
    ascending one: for `_vsOthers` rows the child mean must exceed the parent
    mean; `_2way` rows are kept on the "High" side by construction.
    """
    if splits.empty:
        return pd.DataFrame(columns=["tf", "score", "rank"])
    df = splits.dropna(subset=["score"]).copy()
    up = df.avg_target > df.avg_other
    if direction == "down":
        up = ~up
    if direction == "up":
        up &= (df.kind == "2way") | (df.delta > 0)
    else:
        up &= (df.kind == "2way") | (df.delta < 0)
    df = df[up & (df.n_target_bound.fillna(0) >= min_targets)]
    if df.empty:
        return pd.DataFrame(columns=["tf", "score", "rank"])
    best = df.sort_values("score").groupby("tf", as_index=False).first()
    best["n_branches"] = df.groupby("tf").size().reindex(best.tf).values
    best = best.sort_values("score").reset_index(drop=True)
    best["rank"] = np.arange(1, len(best) + 1)
    return best[["tf", "score", "rank", "n_branches", "coeff", "depth",
                 "parent_mean", "child_mean", "kind", "split"]]
