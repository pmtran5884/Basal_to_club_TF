"""De novo regulon inference with ARACNe-AP over airway metacells."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from basal_to_club.utils.io import load_yaml, setup_logging


def tf_universe(genes, path: str = "data/external/lambert_tfs.txt") -> list[str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"{p} missing. The TF census defines the tested universe and therefore "
            "every rank percentile in the benchmark; it must be explicit, not inferred.")
    tfs = {ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()}
    return sorted({g for g in genes if str(g).upper() in tfs})


def write_aracne_input(mc, path: Path) -> None:
    expr = pd.DataFrame(np.asarray(mc.X.todense() if hasattr(mc.X, "todense") else mc.X),
                        index=mc.obs_names, columns=mc.var_names).T
    expr.index.name = "gene"
    expr.to_csv(path, sep="\t")


def run_aracne(expr_path: Path, tf_path: Path, outdir: Path, p_threshold: float,
               bootstraps: int, seed: int = 0, jar: str | None = None) -> Path:
    jar = jar or shutil.which("aracne3") or "aracne.jar"
    outdir.mkdir(parents=True, exist_ok=True)
    base = ["java", "-Xmx32G", "-jar", jar, "-e", str(expr_path), "-o", str(outdir),
            "--tfs", str(tf_path)]
    subprocess.run(base + ["--pvalue", str(p_threshold), "--seed", str(seed),
                           "--calculateThreshold"], check=True)
    for i in range(bootstraps):
        subprocess.run(base + ["--pvalue", str(p_threshold), "--seed", str(i + 1)], check=True)
    subprocess.run(["java", "-Xmx32G", "-jar", jar, "-o", str(outdir), "--consolidate"],
                   check=True)
    return outdir / "network.txt"


def to_regulon(network: pd.DataFrame, expr: pd.DataFrame, min_targets: int,
               max_targets: int) -> pd.DataFrame:
    """ARACNe gives undirected MI; sign comes from TF-target Spearman, as for the priors."""
    from basal_to_club.network.prior import spearman_signs
    pairs = list(zip(network["regulator"], network["target"], strict=True))
    rho = spearman_signs(expr, pairs)
    out = network.assign(
        mo=np.sign(rho.to_numpy()),
        likelihood=network["mi"] / network.groupby("regulator")["mi"].transform("max"))
    out.loc[out["mo"] == 0, "mo"] = 1.0
    out = out.rename(columns={"regulator": "tf"})
    sizes = out.groupby("tf")["target"].size()
    out = out[out["tf"].isin(sizes[sizes >= min_targets].index)]
    return (out.sort_values(["tf", "likelihood"], ascending=[True, False])
               .groupby("tf", group_keys=False).head(max_targets)
               [["tf", "target", "mo", "likelihood"]].reset_index(drop=True))


def main(sm):
    log = setup_logging(sm.log[0])
    params = load_yaml("config/params.yaml")
    cfg = params["aracne"]
    mc = sc.read_h5ad(sm.input.h5ad)
    if "per_dataset" in sm.wildcards.network:
        log.info("per-dataset ARACNe: %d datasets", mc.obs["dataset_id"].nunique())

    detected = (np.asarray((mc.X > 0).mean(axis=0)).ravel() >= cfg["min_metacell_detection_frac"])
    mc = mc[:, detected].copy()
    expr = pd.DataFrame(np.asarray(mc.X.todense() if hasattr(mc.X, "todense") else mc.X),
                        index=mc.obs_names, columns=mc.var_names)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        write_aracne_input(mc, tmp / "expr.tsv")
        (tmp / "tfs.txt").write_text("\n".join(tf_universe(mc.var_names)))
        net_path = run_aracne(tmp / "expr.tsv", tmp / "tfs.txt", tmp / "out",
                              cfg["p_threshold"], cfg["bootstraps"], params.get("seed", 0))
        network = pd.read_csv(net_path, sep="\t")
    regulon = to_regulon(network, expr, params["regulon_filter"]["min_targets"],
                         params["regulon_filter"]["max_targets"])
    log.info("%s: %d TFs, %d edges", sm.wildcards.network, regulon["tf"].nunique(), len(regulon))
    regulon.to_csv(sm.output[0], sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
