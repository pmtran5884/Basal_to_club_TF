"""Config loading, provenance, and the Snakemake script-mode entry helper."""
from __future__ import annotations

import hashlib
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

LOG = logging.getLogger("basal_to_club")


def setup_logging(logfile: str | None = None, level: int = logging.INFO) -> logging.Logger:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(logfile))
    logging.basicConfig(
        level=level, force=True, handlers=handlers,
        format="%(asctime)s %(levelname)s %(name)s | %(message)s",
    )
    return LOG


def load_yaml(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def load_configs(root: str | Path = ".") -> dict[str, dict]:
    """Load all five configs. `frozen.yaml` overrides `params.yaml` when present."""
    root = Path(root)
    cfg = {name: load_yaml(root / "config" / f"{name}.yaml")
           for name in ("datasets", "markers", "params", "benchmark", "regulons")}
    frozen = root / "config" / "frozen.yaml"
    if frozen.exists():
        cfg["params"] = deep_update(cfg["params"], load_yaml(frozen))
        cfg["params"]["_frozen"] = True
    return cfg


def deep_update(base: dict, overlay: dict) -> dict:
    out = dict(base)
    for k, v in overlay.items():
        out[k] = deep_update(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def dataset_spec(cfg: dict, dataset_id: str) -> dict:
    for d in cfg["datasets"]["datasets"]:
        if d["id"] == dataset_id:
            return d
    raise KeyError(f"{dataset_id} not in config/datasets.yaml")


def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def write_provenance(path: str | Path, record: dict[str, Any]) -> None:
    """Every external download gets one of these. Reproducibility is not optional."""
    record = dict(record)
    record["retrieved_utc"] = datetime.now(timezone.utc).isoformat()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)


def write_json(path: str | Path, obj: Any) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=str)


def snakemake_context():
    """Return the injected `snakemake` object in script mode, else None."""
    import builtins
    return getattr(builtins, "snakemake", None) or globals().get("snakemake")
