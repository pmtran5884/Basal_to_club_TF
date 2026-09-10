"""CollecTRI: curated, already-signed. Used both as a regulon and as the audit
reference for the sign heuristic applied to the ChEA3 libraries."""
from __future__ import annotations

from basal_to_club.utils.io import load_yaml, setup_logging


def main(sm):
    import decoupler as dc
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    net = dc.op.collectri(organism="human") if hasattr(dc, "op") else dc.get_collectri(organism="human")
    net = net.rename(columns={"source": "tf", "weight": "mo"})
    net["mo"] = net["mo"].clip(-1, 1)
    net["likelihood"] = 1.0
    sizes = net.groupby("tf")["target"].size()
    net = net[net["tf"].isin(sizes[sizes >= params["regulon_filter"]["min_targets"]].index)]
    log.info("collectri: %d TFs, %d edges", net["tf"].nunique(), len(net))
    net[["tf", "target", "mo", "likelihood"]].to_csv(sm.output[0], sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
