"""Stage 3: re-parse every DREM model directory into one ranking table.

Separate from the driver so rankings can be rebuilt without refitting models.
Path (edge) tables are preferred when present -- they test enrichment of a TF's
targets among the genes taking an edge -- and the split tables, which DREM
writes from its per-split logistic regression, are the fallback.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from basal_to_club.drem import parse as DP  # noqa: E402

# Which time course to model: "drem_ali" (nasal, GSE121600) or "drem_bronchial"
# (GSE233145). Both prep scripts write ali.h5ad + cell_classes.csv into their
# own results dir, so the grid code is identical for the two datasets.
OUT = ROOT / "results" / os.environ.get("BTC_DREM_SET", "drem_ali")


def main():
    rows, runs = [], []
    for tfdir in sorted((OUT / "models").glob("*/*/*/*/tfscores")):
        arm, donor, lineage, family = tfdir.parts[-5:-1]
        edges = DP.parse_edges(tfdir)
        if not edges.empty:
            rk, source = DP.rank_tfs(edges), "path_tables"
        else:
            rk, source = DP.rank_tfs_from_splits(DP.parse_splits(tfdir)), "split_tables"
        done = (tfdir.parent / "drem_stdout.txt").exists()
        clean = done and "Exception" not in (tfdir.parent / "drem_stdout.txt").read_text()
        runs.append({"arm": arm, "donor": donor, "lineage": lineage, "family": family,
                     "finished": done, "clean_exit": clean, "source": source,
                     "n_path_tables": len(list(tfdir.glob("path_*.txt"))),
                     "n_split_tables": len(list(tfdir.glob("split_*.txt"))),
                     "n_tfs": int(0 if rk is None or rk.empty else rk.tf.nunique())})
        if rk is not None and not rk.empty:
            rk.to_csv(tfdir.parent / "tf_ranking.csv", index=False)
            rows.append(rk.assign(arm=arm, donor=donor, lineage=lineage, family=family))
    runs = pd.DataFrame(runs)
    runs.to_csv(OUT / "drem_runs.csv", index=False)
    if rows:
        allrk = pd.concat(rows, ignore_index=True)
        allrk.to_csv(OUT / "drem_tf_rankings.csv", index=False)
        print(f"{len(allrk)} ranking rows over {allrk.tf.nunique()} TFs")
    json.dump({"n_runs": int(len(runs)), "n_finished": int(runs.finished.sum()),
               "n_clean": int(runs.clean_exit.sum())},
              open(OUT / "drem_collect_summary.json", "w"), indent=1)
    print(runs.to_string(index=False))


if __name__ == "__main__":
    main()
