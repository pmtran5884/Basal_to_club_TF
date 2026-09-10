"""Download raw matrices. Fetch is isolated from parsing so the rest of the DAG
runs offline from cache, and so a re-run never silently re-downloads."""
from __future__ import annotations

from pathlib import Path

from basal_to_club.utils.io import setup_logging, sha256, write_provenance

GEO_SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{acc}/suppl/"
CXG_ASSET = "https://api.cellxgene.cziscience.com/curation/v1/datasets/{ds}/assets"


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


def fetch_cellxgene(dataset_id: str, dest: Path) -> list[Path]:
    """Curated .h5ad with ontology-standardized labels, preferred over reprocessing raw."""
    import requests
    dest.mkdir(parents=True, exist_ok=True)
    assets = requests.get(CXG_ASSET.format(ds=dataset_id), timeout=300).json()
    h5ad = next(a for a in assets if a["filetype"].upper() == "H5AD")
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
        files = fetch_cellxgene(spec["cellxgene_dataset_id"], dest)
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
