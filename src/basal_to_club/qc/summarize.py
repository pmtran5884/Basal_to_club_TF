"""Collate per-dataset QC metrics into one auditable table."""
from __future__ import annotations

import json

import pandas as pd


def main(sm):
    rows = [json.load(open(p)) for p in sm.input]
    df = pd.DataFrame(rows).sort_values("dataset_id")
    df.to_csv(sm.output[0], index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
