"""Concatenate QC'd datasets on the intersection of genes."""
from __future__ import annotations

import anndata as ad

from basal_to_club.utils.io import setup_logging


def main(sm):
    log = setup_logging(sm.log[0])
    parts = [ad.read_h5ad(p) for p in sm.input]
    atlas = ad.concat(parts, join="inner", label="batch_source",
                      keys=[p.obs["dataset_id"].iloc[0] for p in parts],
                      index_unique="-", merge="unique")
    log.info("concatenated %d datasets -> %d cells x %d genes (intersection)",
             len(parts), atlas.n_obs, atlas.n_vars)
    if atlas.n_vars < 8000:
        raise ValueError(
            f"gene intersection collapsed to {atlas.n_vars}. Almost always an "
            "identifier mismatch (Ensembl IDs vs symbols, or unresolved aliases) "
            "rather than real biology.")
    atlas.write_h5ad(sm.output[0], compression="gzip")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
