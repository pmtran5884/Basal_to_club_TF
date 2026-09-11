"""Headless DREM execution."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def java_binary() -> str:
    """The conda JDK, which is not on PATH as `java` on macOS (the system stub is)."""
    prefix = os.environ.get("CONDA_PREFIX", "")
    cand = Path(prefix) / "lib/jvm/bin/java"
    if cand.exists():
        return str(cand)
    found = shutil.which("java")
    if found is None:
        raise RuntimeError("no Java runtime found; install openjdk into the environment")
    return found


def run_drem(jar: str | Path, settings: str | Path, out_model: str | Path,
             gene_assign: str | Path, tf_dir: str | Path,
             workdir: str | Path = ".", heap: str = "4g",
             timeout: int = 3600,
             classes_dir: str | Path | None = None) -> subprocess.CompletedProcess:
    """Run one DREM batch job. Returns the completed process (stdout captured).

    DREM lists a `TFInput` directory relative to the working directory at
    startup; it only warns when absent, but creating it keeps the log clean.
    """
    workdir = Path(workdir)
    (workdir / "TFInput").mkdir(parents=True, exist_ok=True)
    Path(tf_dir).mkdir(parents=True, exist_ok=True)
    jar = Path(jar).resolve()
    if classes_dir is None:
        launch = ["-jar", str(jar)]
    else:
        # DREM_IO_Batch is recompiled locally with a null-guard on the child
        # ordering array (see scripts/patch_drem.py); it must precede the jar.
        launch = ["-cp", f"{Path(classes_dir).resolve()}:{jar}",
                  "edu.cmu.cs.sb.drem.DREM_IO"]
    cmd = [java_binary(), f"-Xmx{heap}", "-Djava.awt.headless=true", *launch,
           "-b", str(Path(settings).resolve()),
           str(Path(out_model).resolve()), str(Path(gene_assign).resolve()),
           str(Path(tf_dir).resolve())]
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=timeout)
    return proc


def drem_succeeded(proc: subprocess.CompletedProcess, tf_dir: str | Path) -> bool:
    """DREM exits 0 even on a caught exception, so check for actual output."""
    if proc.returncode != 0:
        return False
    if "Exception" in (proc.stdout or "") or "Error" in (proc.stdout or ""):
        return False
    return any(Path(tf_dir).glob("path_*.txt"))
