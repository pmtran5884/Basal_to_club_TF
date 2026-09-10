rule fetch_dataset:
    """Download raw matrices. Isolated so the rest of the DAG can run offline from cache."""
    output:
        raw=protected("data/raw/{dataset}/.fetched"),
    params:
        spec=lambda w: next(d for d in DATASETS["datasets"] if d["id"] == w.dataset),
    log:
        "logs/fetch/{dataset}.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/ingest/fetch.py"


rule load_dataset:
    """Raw -> AnnData with harmonized obs schema and a provenance record.

    Sample-level exclusions (non-human samples in GSE121600, infected cells in
    GSE166766) are asserted HERE, not left to a downstream subset.
    """
    input:
        "data/raw/{dataset}/.fetched",
    output:
        h5ad="data/interim/{dataset}_loaded.h5ad",
        prov="results/provenance/{dataset}_load.json",
    params:
        spec=lambda w: next(d for d in DATASETS["datasets"] if d["id"] == w.dataset),
    log:
        "logs/load/{dataset}.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/ingest/load.py"
