"""GSE233145 -- primary human bronchial epithelial ALI differentiation time course.

Wolf et al., serial single-cell sampling of pHBEC cultures at ALI days
0/3/5/7/14/21/28 from two pulmonary-healthy and two end-stage COPD donors.
This is the only human *bronchial*, non-infection, single-cell ALI time course
in GEO with >=3 sampled differentiation days (see
docs/drem_bronchial/DATASET_SCREEN.md).

Two properties of the deposit dictate this loader:

1. The per-sample matrices (`GSM*_<id>.dge.txt.gz`) are Drop-seq digital
   expression tables -- dense tab-separated text, genes x barcodes, *before*
   cell calling (8k-42k barcode columns per sample against ~2.4k published
   cells). QC is therefore ours to redo; it is done per sample with the same
   thresholds as the nasal arm so the two datasets are comparable.
2. `GSE233145_cells_metadata.txt.gz` carries the published per-cell
   annotation (patient, time_point, health_state, cell_type) but *no barcode
   column*, so those labels cannot be joined to the matrices cell by cell.
   They are used as an independent check on composition instead
   (`published_composition`).
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CHUNK_GENES = 2000
# Drop-seq is far shallower than the 10x nasal arm; these are the deposit's own
# thresholds, read off GSE233145_cells_metadata.txt.gz (min n_counts 300,
# min n_genes 200, max percent_mito 0.20). Using the 10x thresholds instead
# discards ~75% of the published cells.
MIN_COUNTS = 300
MIN_GENES = 200
MAX_PCT_MITO = 20.0


def sample_table(path: Path) -> pd.DataFrame:
    """Sample sheet built from the GEO sample titles ('Donor 1, day 14')."""
    return pd.read_csv(path, sep="\t")


def read_dge(path: Path) -> "sc.AnnData":
    """One Drop-seq DGE text matrix -> cells x genes AnnData (raw counts).

    Read in gene blocks and converted to sparse per block: the full dense table
    is ~3e8 values per sample, the sparse form is ~2 orders of magnitude smaller.
    """
    blocks, genes, cells = [], [], None
    dtypes = defaultdict(lambda: np.int32)     # gene column is text, counts are int
    dtypes["GENE"] = str
    for chunk in pd.read_csv(path, sep="\t", index_col=0, chunksize=CHUNK_GENES,
                             dtype=dtypes, compression="gzip"):
        if cells is None:
            cells = list(chunk.columns)
        blocks.append(sparse.csr_matrix(chunk.to_numpy()))
        genes += list(chunk.index)
    X = sparse.vstack(blocks).T.tocsr()          # cells x genes
    ad = sc.AnnData(X, obs=pd.DataFrame(index=cells),
                    var=pd.DataFrame(index=pd.Index(genes, name=None)))
    ad.var_names_make_unique()
    return ad


def qc(ad) -> "sc.AnnData":
    sc.pp.filter_cells(ad, min_counts=MIN_COUNTS)
    sc.pp.filter_cells(ad, min_genes=MIN_GENES)
    mt = [g for g in ad.var_names if g.startswith("MT-")]
    tot = np.asarray(ad.X.sum(1)).ravel()
    pmt = np.asarray(ad[:, mt].X.sum(1)).ravel() / np.maximum(tot, 1) * 100 if mt else np.zeros(ad.n_obs)
    ad.obs["pct_mito"] = pmt
    return ad[pmt < MAX_PCT_MITO].copy()


def load(raw_dir: Path, samples: pd.DataFrame, health_state: str | None = "healthy",
         log=print) -> "sc.AnnData":
    """Concatenated, QC'd, log-normalised object with `counts` layer.

    obs: sample_id, gsm, donor, day, health_state. Replicate samples of the same
    donor-day keep their own sample_id and are pooled downstream by day.
    """
    sel = samples if health_state is None else samples[samples.health_state == health_state]
    ads = []
    for r in sel.itertuples():
        a = qc(read_dge(Path(raw_dir) / r.file))
        a.obs["sample_id"] = r.sample_id
        a.obs["gsm"] = r.gsm
        a.obs["donor"] = r.donor
        a.obs["day"] = int(r.day)
        a.obs["health_state"] = r.health_state
        a.obs["sample"] = f"{r.donor}_d{r.day}" + ("" if r.replicate == 1 else f"_r{r.replicate}")
        a.obs_names = [f"{r.sample_id}_{b}" for b in a.obs_names]
        ads.append(a)
        log(f"{r.sample_id} ({r.donor} day {r.day}): {a.n_obs} cells kept of {r.n_barcodes_raw} barcodes")
    # outer, not inner: a DGE file lists only the genes detected in that sample,
    # so "absent" means zero counts, not unmeasured. An inner join silently drops
    # any gene undetected in even one sample -- it removed FOXJ1 and MUC5B (absent
    # from the day-0 samples, which have no ciliated or mucous cells), i.e. exactly
    # the markers and control TFs the analysis is built on.
    ad = sc.concat(ads, join="outer", fill_value=0)
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad


def published_composition(metadata_gz: Path) -> pd.DataFrame:
    """Per-sample cell-type fractions as annotated by the depositing authors.

    The published classes are Basal_1, Basal_2, Suprabasal, Secretory,
    Transitional Ciliated and Ciliated -- note there is no club/goblet split,
    which is what the resolver has to supply.
    """
    md = pd.read_csv(metadata_gz, sep="\t", compression="gzip")
    tab = md.groupby(["sample_id", "cell_type"]).size().unstack(fill_value=0)
    return tab.div(tab.sum(1), axis=0)
