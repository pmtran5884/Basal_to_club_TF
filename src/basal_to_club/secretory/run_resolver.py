"""Apply the club/goblet resolver to the secretory compartment."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from basal_to_club.secretory.resolver import resolve_secretory
from basal_to_club.utils.io import load_yaml, setup_logging, write_json

SECRETORY_LABELS = ("secretory", "club", "goblet", "club/goblet", "secretory_unassigned",
                    "mucous", "serous", "hybrid")


def main(sm):
    log = setup_logging(sm.log[0])
    params, markers = load_yaml("config/params.yaml"), load_yaml(sm.input.markers)
    adata = sc.read_h5ad(sm.input.h5ad)

    coarse = adata.obs["cell_type_coarse"].astype(str).str.lower()
    is_sec = coarse.str.contains("|".join(SECRETORY_LABELS))
    log.info("secretory compartment: %d/%d cells enter the resolver", is_sec.sum(), adata.n_obs)
    if is_sec.sum() < 500:
        raise ValueError(
            f"only {is_sec.sum()} secretory cells found. Either the coarse annotation "
            "vocabulary is unexpected or the epithelial gate is too aggressive; "
            "resolving a compartment this small would not be meaningful.")

    sub = adata[is_sec].copy()
    assignments, thresholds = resolve_secretory(sub, markers["markers"], params)

    adata.obs["cell_class"] = adata.obs["cell_type_coarse"].astype(str)
    adata.obs.loc[assignments.index, "cell_class"] = assignments["cell_class"].astype(str)
    for col in ("club", "goblet", "ciliated", "smg"):
        adata.obs[f"resolver_score_{col}"] = assignments[col].reindex(adata.obs_names)
    adata.obs["cell_class"] = adata.obs["cell_class"].astype("category")
    adata.uns["resolver_thresholds"] = thresholds

    counts = adata.obs["cell_class"].value_counts().to_dict()
    log.info("resolved classes: %s", counts)
    write_json("results/secretory/class_counts.json", {k: int(v) for k, v in counts.items()})
    assignments.to_csv(sm.output.assignments)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(assignments["club"], assignments["goblet"], s=2, alpha=0.3,
               c=assignments["cell_class"].cat.codes, cmap="tab10")
    ax.set_xlabel("club program score")
    ax.set_ylabel("goblet program score")
    ax.set_title("Club/goblet separation after SMG removal")
    fig.savefig(sm.output.fig, dpi=200, bbox_inches="tight")
    adata.write_h5ad(sm.output.h5ad, compression="gzip")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
