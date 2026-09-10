"""Download raw matrices. Fetch is isolated from parsing so the rest of the DAG
runs offline from cache, and so a re-run never silently re-downloads."""
from __future__ import annotations

from pathlib import Path

from basal_to_club.utils.io import setup_logging, sha256, write_provenance

GEO_SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{acc}/suppl/"
CXG_ASSET = "https://api.cellxgene.cziscience.com/curation/v1/datasets/{ds}/assets"
CXG_COLLECTION = "https://api.cellxgene.cziscience.com/curation/v1/collections/{col}"


def fetch_geo(accession: str, dest: Path) -> list[Path]:
    import requests
    dest.mkdir(parents=True, exist_ok=True)
    url = GEO_SUPP.format(stub=accession[:-3], acc=accession)
    listing = requests.get(url, timeout=300)
    listing.raise_for_status()
    names = sorted(set(__import__("re").findall(r'href="([^"]+\.(?:tar|gz|mtx|h5|h5ad|csv|txt))"',
                                                listing.text)))
    out = []
    for name in names:
        target = dest / name
        if target.exists():
            out.append(target)
            continue
        with requests.get(url + name, stream=True, timeout=1800) as r:
            r.raise_for_status()
            with open(target, "wb") as fh:
                for chunk in r.iter_content(1 << 22):
                    fh.write(chunk)
        out.append(target)
    return out


def _cxg_h5ad_asset(dataset_id: str, collection_id: str | None) -> dict:
    """Resolve the H5AD asset for a dataset.

    The bare /datasets/{id}/assets endpoint returns 404 as of 2026-09; the
    collection-scoped listing is the supported path. The bare endpoint is still
    tried first so that a future restoration of it keeps working, and its error
    payload (a dict, not a list) is detected rather than being indexed blindly.
    """
    import requests
    bare = requests.get(CXG_ASSET.format(ds=dataset_id), timeout=300)
    if bare.ok:
        payload = bare.json()
        if isinstance(payload, list):
            return next(a for a in payload if a["filetype"].upper() == "H5AD")
    if not collection_id:
        raise RuntimeError(
            f"CELLxGENE dataset {dataset_id}: the per-dataset asset endpoint returned "
            f"{bare.status_code} and no cellxgene_collection_id is configured for this "
            "dataset, so the asset URL cannot be resolved. Add "
            "`cellxgene_collection_id` to its entry in config/datasets.yaml."
        )
    col = requests.get(CXG_COLLECTION.format(col=collection_id), timeout=300)
    col.raise_for_status()
    hits = [d for d in col.json().get("datasets", []) if d["dataset_id"] == dataset_id]
    if len(hits) != 1:
        raise RuntimeError(
            f"dataset {dataset_id} matched {len(hits)} entries in collection "
            f"{collection_id}; expected exactly 1"
        )
    return next(a for a in hits[0]["assets"] if a["filetype"].upper() == "H5AD")


def fetch_cellxgene(dataset_id: str, dest: Path,
                    collection_id: str | None = None) -> list[Path]:
    """Curated .h5ad with ontology-standardized labels, preferred over reprocessing raw."""
    import requests
    dest.mkdir(parents=True, exist_ok=True)
    h5ad = _cxg_h5ad_asset(dataset_id, collection_id)
    target = dest / f"{dataset_id}.h5ad"
    if not target.exists():
        with requests.get(h5ad["url"], stream=True, timeout=7200) as r:
            r.raise_for_status()
            with open(target, "wb") as fh:
                for chunk in r.iter_content(1 << 22):
                    fh.write(chunk)
    return [target]


def main(sm):
    log = setup_logging(sm.log[0])
    spec = sm.params.spec
    dest = Path(f"data/raw/{spec['id']}")
    if spec["source"] == "geo":
        files = fetch_geo(spec["accession"], dest)
    elif spec["source"] == "cellxgene":
        files = fetch_cellxgene(spec["cellxgene_dataset_id"], dest,
                                spec.get("cellxgene_collection_id"))
    else:
        raise ValueError(f"unknown source {spec['source']}")
    write_provenance(f"results/provenance/{spec['id']}_fetch.json", {
        "dataset_id": spec["id"], "source": spec["source"],
        "accession": spec.get("accession"), "citation": spec.get("citation"),
        "files": [{"name": f.name, "sha256": sha256(f), "bytes": f.stat().st_size} for f in files]})
    log.info("fetched %d files for %s", len(files), spec["id"])
    Path(sm.output.raw).touch()


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
