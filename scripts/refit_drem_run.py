"""Refit one DREM model directory in place.

Used when a single model crashed in the batch writer and the fix is a patched
DREM class: refitting the whole grid is hours of compute for one cell.

    python scripts/refit_drem_run.py atlas_label D275 goblet binding
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from basal_to_club.drem import run as DR  # noqa: E402

JAR = ROOT.parent / "external/STEM_DREM/drem.jar"
CLASSES = ROOT.parent / "external/STEM_DREM/patched"


def main(arm, donor, lineage, family):
    d = ROOT / "results/drem_ali/models" / arm / donor / lineage / family
    assert (d / "settings.txt").exists(), f"no prepared inputs in {d}"
    for f in (d / "tfscores").glob("*.txt"):
        f.unlink()
    t0 = time.time()
    proc = DR.run_drem(JAR, d / "settings.txt", d / "model.txt", d / "geneassign.txt",
                       d / "tfscores", workdir=d, heap="4g", timeout=36000,
                       classes_dir=CLASSES)
    (d / "drem_stdout.txt").write_text((proc.stdout or "") + (proc.stderr or ""))
    n_path = len(list((d / "tfscores").glob("path_*.txt")))
    n_split = len(list((d / "tfscores").glob("split_*.txt")))
    print(f"{arm}/{donor}/{lineage}/{family}: {round(time.time()-t0)}s "
          f"path={n_path} split={n_split} "
          f"clean={'Exception' not in (d / 'drem_stdout.txt').read_text()}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
