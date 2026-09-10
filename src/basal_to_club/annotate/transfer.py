"""Label transfer from the HLCA core reference via scArches/scANVI surgery."""
from __future__ import annotations

import numpy as np
import scanpy as sc

from basal_to_club.utils.io import load_yaml, setup_logging


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    query, ref = sc.read_h5ad(sm.input.query), sc.read_h5ad(sm.input.reference)
    try:
        import scvi
        scvi.settings.seed = params.get("seed", 0)
        shared = ref.var_names.intersection(query.var_names)
        ref, query = ref[:, shared].copy(), query[:, shared].copy()
        label_key = next(k for k in ("ann_finest_level", "cell_type", "author_cell_type")
                         if k in ref.obs)
        scvi.model.SCVI.setup_anndata(ref, layer="counts", batch_key="sample_id")
        base = scvi.model.SCVI(ref, n_latent=params["integrate"]["n_latent"])
        base.train(max_epochs=params["integrate"]["max_epochs"])
        scan = scvi.model.SCANVI.from_scvi_model(base, unlabeled_category="unknown",
                                                 labels_key=label_key)
        scan.train(max_epochs=100)
        surgery = scvi.model.SCANVI.load_query_data(query, scan)
        surgery.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})
        query.obs["transferred_label"] = surgery.predict()
        probs = surgery.predict(soft=True)
        query.obs["transfer_confidence"] = np.asarray(probs).max(axis=1)
    except Exception as exc:                            # noqa: BLE001
        # Fail visibly: an unlabelled atlas must not look like a successfully
        # transferred one two rules downstream.
        log.error("label transfer failed: %s", exc)
        query.obs["transferred_label"] = "transfer_failed"
        query.obs["transfer_confidence"] = 0.0
        query.uns["transfer_error"] = str(exc)
    query.write_h5ad(sm.output[0], compression="gzip")


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
