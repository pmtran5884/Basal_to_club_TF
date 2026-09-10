"""Download ChEA3 gene-set libraries once, with provenance.

The web API is deliberately not used at analysis time: a versioned local copy
with a recorded sha256 is what makes a run reproducible six months later.
Confirm the current library names against the ChEA3 site at build time.
"""
from __future__ import annotations

from pathlib import Path

from basal_to_club.utils.io import load_yaml, setup_logging, sha256, write_provenance


def main(sm):
    import requests
    log = setup_logging(sm.log[0])
    cfg = load_yaml("config/regulons.yaml")["prior"]["chea3"]
    dest = Path(sm.output[0])
    dest.mkdir(parents=True, exist_ok=True)

    records = []
    for lib in cfg["libraries"]:
        url = f"{cfg['base_url']}/{lib}.gmt"
        target = dest / f"{lib}.gmt"
        if not target.exists():
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            target.write_bytes(r.content)
        records.append({"library": lib, "url": url, "sha256": sha256(target),
                        "bytes": target.stat().st_size,
                        "n_sets": sum(1 for _ in open(target))})
        log.info("chea3 %s: %d sets", lib, records[-1]["n_sets"])
    write_provenance(sm.output.prov, {"source": "ChEA3", "libraries": records})


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
