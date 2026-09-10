import sys, shutil
from pathlib import Path
sys.path.insert(0, "src")
from basal_to_club.ingest.fetch import fetch_cellxgene
DS = "066943a2-fdac-4b29-b348-40cede398e4e"
dest = Path("data/raw/hlca_core")
[p] = fetch_cellxgene(DS, dest, "6f6d381a-7701-4781-935c-db10d30de293")
final = dest / "hlca_core.h5ad"
if not final.exists():
    p.rename(final)
print("ready", final, final.stat().st_size / 1e9, "GB", flush=True)
